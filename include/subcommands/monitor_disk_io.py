#!/usr/bin/env python

# Run this tool in two stages:
#
# 1. Collect data: ./monitor-disk-io.py [command]
# 2. Print statistics: ./monitor-disk-io.py
# 3. $$$
#
# The first command will run strace and write the output to
# _monitor_disk_io/strace-out.txt. The second call (without any arguments)
# will parse that file and try to find out which file every read, write,
# and seek corresponds to.
# If you want to profile a BASH command, put it in a script and pass the
# path of that script (i. e. BASH voodoo is not allowed here but can be
# accomplished if a BASH script is used).

import yaml
import re
import os
import subprocess
import logging
import glob
import copy
import sys
WRITE_PROC_FILES = False


logger = logging.getLogger("uap_logger")

proc_files = {}

if not os.path.exists('_monitor_disk_io'):
    os.mkdir('_monitor_disk_io')

if len(sys.argv) > 1:
    pigz = subprocess.Popen(
        "pigz -p 2 -b 4096 -c > _monitor_disk_io/strace-out.txt.gz",
        stdin=subprocess.PIPE,
        shell=True)
    args = ["strace", "-f", "-o", '/dev/stderr']
    args.extend(sys.argv[1:])
    p = subprocess.Popen(args, stderr=subprocess.PIPE)
    strace_out = p.stderr
    for line in strace_out:
        pigz.stdin.write(line)
    pigz.stdin.close()
    pigz.wait()
    exit(0)

pigz = subprocess.Popen("pigz -p 1 -d -c _monitor_disk_io/strace-out.txt.gz",
                        stdout=subprocess.PIPE, shell=True)
strace_out = pigz.stdout

if len(glob.glob('_monitor_disk_io/*.proc.txt')) > 0:
    os.system("rm _monitor_disk_io/*.proc.txt")

path_for_pid_and_fd = {}
stats = {}


def handle_line(pid, line):
    line = line.strip()
    if WRITE_PROC_FILES:
        if pid not in proc_files:
            proc_files[pid] = open("_monitor_disk_io/%s.proc.txt" % pid, 'w')
        proc_files[pid].write(line + "\n")
    m = re.search(r'^(\w+)\((.*)\)\s+=\s+(.+)$', line)
    if m:
        command = str(m.group(1))
        args = str(m.group(2)).strip()
        retval = str(m.group(3))
        if command == 'clone':
            for _ in path_for_pid_and_fd[pid].keys():
                if retval not in path_for_pid_and_fd:
                    path_for_pid_and_fd[retval] = {}
                path_for_pid_and_fd[retval][_] = copy.copy(
                    path_for_pid_and_fd[pid][_])

        if command == 'dup2':
            fds = [_.strip() for _ in args.split(',')]
            try:
                path_for_pid_and_fd[pid][fds[1]] = copy.copy(
                    path_for_pid_and_fd[pid][fds[0]])
            except BaseException:
                path_for_pid_and_fd[pid][fds[1]] = '[unknown]'

        if command == 'open':
            if pid not in path_for_pid_and_fd:
                path_for_pid_and_fd[pid] = {}
            path_for_pid_and_fd[pid][retval] = re.search(
                "^\\\"([^\\\"]+)\\\"", args).group(1)
        if command == 'close':
            if pid not in path_for_pid_and_fd:
                path_for_pid_and_fd[pid] = {}
            fd = args.strip()
        if command == 'lseek':
            fd = None
            m = re.search(r"^(\d+),", args)
            if m:
                fd = m.group(1)
            if fd:
                path = '[unknown]'
                try:
                    path = path_for_pid_and_fd[pid][fd]
                except BaseException:
                    if fd == '0':
                        path = 'stdin'
                    elif fd == '1':
                        path = 'stdout'
                    elif fd == '2':
                        path = 'stderr'
                    else:
                        pass
                if path not in stats:
                    stats[path] = {'read': {}, 'write': {}, 'lseek': 0}
                stats[path]['lseek'] += 1
        if command == 'read' or command == 'write':
            fd = None
            size = None
            m = re.search(r"^(\d+),", args)
            if m:
                fd = m.group(1)
            size = retval
            if fd and size:
                try:
                    sizek = int(size) / 1024
                except ValueError:
                    return
                path = '[unknown]'
                try:
                    path = path_for_pid_and_fd[pid][fd]
                except BaseException:
                    if fd == '0':
                        path = 'stdin'
                    elif fd == '1':
                        path = 'stdout'
                    elif fd == '2':
                        path = 'stderr'
                    else:
                        pass
                if path not in stats:
                    stats[path] = {'read': {}, 'write': {}, 'lseek': 0}
                if sizek not in stats[path][command]:
                    stats[path][command][sizek] = 0
                stats[path][command][sizek] += 1


def size_to_cat(s):
    if s < 32:
        return (0, '< 32k ')
    elif s < 128:
        return (1, '32k+ ')
    elif s < 1024:
        return (2, '128k+ ')
    elif s < 2048:
        return (3, '1024k+ ')
    elif s < 3072:
        return (4, '2048k+ ')
    elif s < 4096:
        return (5, '3072k+ ')
    elif s < 8192:
        return (6, '4096k+ ')
    else:
        return (7, '8192k+')


line_buffer = {}
for line in strace_out:
    line = line.strip()
    pid = str(re.search(r'^(\d+)\s', line).group(1))
    line = line[line.index(' ') + 1:]
    if 'resumed>' in line:
        line = re.sub(r'\<.+\>', '', line)
        line_buffer[pid] += line
        handle_line(pid, line_buffer[pid])
    else:
        if '<unfinished' in line:
            line = re.sub(r'\<.+\>', '', line)
            line_buffer[pid] = line
        else:
            line_buffer[pid] = line
            handle_line(pid, line_buffer[pid])

for path in stats.keys():
    cancel = False
    for _ in [
        'python_env',
        '/proc',
        '/etc',
        '/usr',
        '.git',
        '.so',
        '.py',
        '.pyc',
        'stdin',
        'stdout',
        'stderr',
            '/dev']:
        if _ in path:
            cancel = True
            continue
    if cancel:
        continue
    printed_path = False
    for mode in ['read', 'write']:
        printed_mode = False
        hist = {}
        mod_size = {}
        for _, count in stats[path][mode].items():
            cat = size_to_cat(_)
            if cat not in mod_size:
                mod_size[cat] = 0
            mod_size[cat] += count
        for key in sorted(mod_size.keys(), reverse=True):
            size = key[1]
            # if key[0] == 0:
            # continue
            if not printed_path:
                print('-' * len(path))
                print(path)
                print('-' * len(path))
                printed_path = True
            if not printed_mode:
                print(mode.upper() + 'S:')
                printed_mode = True
            print('{:>8} {:>5}x'.format(str(size), str(mod_size[key])))
    if stats[path]['lseek'] > 0:
        if not printed_path:
            print('-' * len(path))
            print(path)
            print('-' * len(path))
            printed_path = True
        print("LSEEKS:   " + str(stats[path]['lseek']) + 'x')

for pid, f in proc_files.items():
    f.close()
