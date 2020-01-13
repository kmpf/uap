'''
This module can be used to launch child processes and wait for them.
Processes may either run on their own or pipelines can be built with them.
'''

import sys
sys.path.append('./include/steps')
import copy
import datetime
import errno
import fcntl
import hashlib
from logging import getLogger
import misc
import os
import psutil
import signal
import subprocess
import tempfile
import time
import traceback
import yaml

logger=getLogger("uap_logger")

class TimeoutException(Exception):
    pass

def timeout_handler(signum, frame):
    time.sleep(3600)
    raise TimeoutException()

def restore_sigpipe_handler():
    # http://www.chiark.greenend.org.uk/ucgi/~cjwatson/blosxom/2009-07-02-python-sigpipe.html
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    os.setsid()

class ProcessPool(object):
    '''
    The process pool provides an environment for launching and monitoring
    processes. You can launch any number of unrelated processes plus any number
    of pipelines in which several processes are chained together.

    Use it like this::

        with process_pool.ProcessPool(self) as pool:
            # launch processes or create pipelines here

    When the scope opened by the *with* statement is left, all processes are
    launched and being watched. The process pool then waits until all processes
    have finished. You cannot launch a process pool within another process pool,
    but you can launch multiple pipeline and independent processes within a
    single process pool. Also, you can launch several process pools sequentially.
    '''

    TAIL_LENGTH = 1024
    '''
    Size of the tail which gets recorded from both *stdout* and *stderr* streams
    of every process launched with this class, in bytes.
    '''

    COPY_BLOCK_SIZE = 4194304
    '''
    When *stdout* or *stderr* streams should be written to output files, this is
    the buffer size which is used for writing.
    '''

    SIGTERM_TIMEOUT = 10
    '''
    After a SIGTERM signal is issued, wait this many seconds before going postal.
    '''

    process_watcher_pid = None

    current_instance = None
    process_pool_is_dead = False

    # signal names for numbers... kudos to
    # http://stackoverflow.com/questions/2549939/get-signal-names-from-numbers-in-python
    SIGNAL_NAMES = dict((getattr(signal, n), n) for n in dir(signal) if n.startswith('SIG') and '_' not in n )

    class Pipeline(object):
        '''
        This class can be used to chain multiple processes together.

        Use it like this::

            with pool.Pipeline(pool) as pipeline:
                # append processes to the pipeline here
        '''
        def __init__(self, pool):
            pool.launch_calls.append(self)
            self.append_calls = []

        def __enter__(self):
            return self

        def __exit__(self, type, value, traceback):
            pass

        def append(self, args, stdout_path = None, stderr_path = None, hints = {}):
            '''
            Append a process to the pipeline. Parameters get stored and are passed
            to *ProcessPool.launch()* later, so the same behaviour applies.
            '''
            call = {
                'args': copy.deepcopy(args),
                'stdout_path': copy.copy(stdout_path),
                'stderr_path': copy.copy(stderr_path),
                'hints': copy.deepcopy(hints)
            }
            self.append_calls.append(call)

    def __init__(self, run):
        if ProcessPool.process_pool_is_dead:
            logger.error("We have encountered an error, stopping now...")
            sys.exit(1)
        # the run for which this ProcessPool computes stuff
        # (for temporary paths etc.)
        self._run = run

        # log entries
        self.log_entries = []

        # dict of PID -> path to copy process report
        self.copy_process_reports = {}

        # set of currently running PIDs
        # whenever a child process exits, its PID will get removed from this
        self.running_procs = set()

        # list of PIDs in the order that processes were launched
        # whenever a child process exits, its PID will remain in here
        self.proc_order = []

        # we must keep Popen objects around, or their destructor will be called
        self.popen_procs = {}

        # dict of PID -> process info
        # whenever a child process exits, its PID will remain in here
        self.proc_details = {}

        self.process_watcher_report = dict()

        # list of temp paths to clean up
        self.temp_paths = []

        # list of commands to be launched
        self.launch_calls = []

        # List of processes we killed deliberately. Look: every time a
        # within a pipeline exits, we SIGTERM its predecessor. This is
        # necessary because otherwise, stuff is hanging forever.
        self.ok_to_fail = set()

        self.copy_processes_for_pid = dict()

        self.clean_up = False

    def clean_up_temp_paths(self):
        self.clean_up = True

    def get_run(self):
        return self._run

    def check_subprocess_command(self, command):
        for argument in command:
            if not isinstance(argument, str):
                logger.error(
                    "The command to be launched '%s' " % command +
                    "contains non-string argument '%s'. " % argument +
                    "Therefore the command will fail. Please " +
                    "fix this type issue.")
                sys.exit(1)
        return

    def load_unload_module(self, module_cmd):
        if module_cmd.__class__ == str:
            module_cmd = [module_cmd]

        for command in module_cmd:
            if type(command) is str:
                command = command.split()
            self.check_subprocess_command(command)

            try:
                proc = subprocess.Popen(
                    command,
                    stdin = None,
                    stdout = subprocess.PIPE,
                    stderr = subprocess.PIPE,
                    close_fds = True)

            except OSError as e:
                logger.error("Error while executing '%s' "
                             "Error no.: %s Error message: %s" %
                             (" ".join(command), e.errno, e.strerror))
                sys.exit(1)

            (output, error) = proc.communicate()
            exec output
            sys.stderr.write(error)
            sys.stderr.flush()

        return

    def launch_pre_post_command(self, commands):
        if commands.__class__ == str:
            commands = [commands]

        for command in commands:
            if type(command) is str:
                command = command.split()

            self.launch(command)

    def __enter__(self):
        if ProcessPool.current_instance is not None:
            logger.error("Sorry, only one instance of ProcessPool allowed at "
                         "a time.")
            sys.exit(1)
        ProcessPool.current_instance = self

        # First we have to add the pre_command commands for execution
        pre_commands = self.get_run().get_step().get_pre_commands().values()
        if len(pre_commands) > 0:
            self.launch_pre_post_command(pre_commands)

        return self

    def __exit__(self, type, value, traceback):
        # Lastly we have to add the post_command commands for execution
        post_commands = self.get_run().get_step().get_post_commands().values()
        if len(post_commands) > 0:
            self.launch_pre_post_command(post_commands)

        # before everything is launched load the necessary modules
        module_loads = self.get_run().get_step().get_module_loads().values()
        if len(module_loads) > 0:
            for module_load in module_loads:
                self.load_unload_module(module_load)

        # now launch all processes...
        self._launch_all_processes()

        # ...and wait until all child processes have exited
        try:
            self._wait()
        except:
            # pass log to step even if there was a problem
            self.get_run().get_step().append_pipeline_log(self.get_log())
            raise

        # if there was no exception, still pass log to step
        self.get_run().get_step().append_pipeline_log(self.get_log())

        # after finishing gracefully unload modules
        module_unloads = self.get_run().get_step().get_module_unloads().values()
        if len(module_unloads) > 0:
            for module_unload in module_unloads:
                self.load_unload_module(module_unload)

        # remove all temporary files or directories we know of
        if self.clean_up:
            self.get_run().remove_temporary_paths()

        ProcessPool.current_instance = None

    def launch(self, args, stdout_path = None, stderr_path = None, hints = {},
            working_directory = None):
        '''
        Launch a process. Arguments, including the program itself, are passed in
        *args*. If the program is not a binary but a script which cannot be
        invoked directly from the command line, the first element of *args* must
        be a list like this: *['python', 'script.py']*.

        Use *stdout_path* and *stderr_path* to redirect *stdout* and *stderr*
        streams to files. In any case, the output of both streams gets watched,
        the process pool calculates SHA1 checksums automatically and also keeps
        the last 1024 bytes of every stream. This may be useful if a process
        crashes and writes error messages to *stderr* in which case you can see
        them even if you didn't redirect *stderr* to a log file.

        Hints can be specified but are not essential. They help to determine the
        direction of arrows for the run annotation graphs rendered by GraphViz
        (sometimes, it's not clear from the command line whether a certain file
        is an input or output file to a given process).
        '''
        call = {
            'args': copy.deepcopy(args),
            'stdout_path': copy.copy(stdout_path),
            'stderr_path': copy.copy(stderr_path),
            'working_directory': copy.copy(working_directory),
            'hints': copy.deepcopy(hints)
        }

        self.launch_calls.append(call)

    def log(self, message):
        '''
        Append a message to the pipeline log.
        '''
        formatted_message = "[%s] %s" % (datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), message)
        #sys.stderr.write(formatted_message + "\n")
        self.log_entries.append(formatted_message)

    def get_log(self):
        '''
        Return the log as a dictionary.
        '''
        log = dict()
        log['processes'] = []
        proc_details_copy = dict()
        for pid in self.proc_order:
            if 'listener_for' in self.proc_details[pid]:
                attach_info = self.proc_details[pid]['listener_for']
                proc_details_copy[attach_info[0]][attach_info[1] + '_copy'] = copy.deepcopy(self.proc_details[pid])
                del proc_details_copy[attach_info[0]][attach_info[1] + '_copy']['listener_for']
            else:
                proc_details_copy[pid] = copy.deepcopy(self.proc_details[pid])

        for pid in self.proc_order:
            if pid in proc_details_copy:
                log['processes'].append(proc_details_copy[pid])

        log['log'] = copy.deepcopy(self.log_entries)
        log['process_watcher'] = copy.deepcopy(self.process_watcher_report)
        log['ok_to_fail'] = copy.deepcopy(self.ok_to_fail)

        return log

    def _launch_all_processes(self):
        for info in self.launch_calls:
            if info.__class__ == ProcessPool.Pipeline:
                pipeline = info
                use_stdin = None
                last_pid = None
                for index, info in enumerate(pipeline.append_calls):
                    use_stdin, pid = self._do_launch(info,
                        index < len(pipeline.append_calls) - 1, use_stdin)
                    if last_pid is not None:
                        self.proc_details[pid]['use_stdin_of'] = last_pid
                    last_pid = pid
            else:
                self._do_launch(info)

    def _do_launch(self, info, keep_stdout_open = False, use_stdin = None):
        '''
        Launch a process and after that, launch a copy process for *stdout* and
        *stderr* each.
        '''
        args = copy.deepcopy(info['args'])
        stdout_path = copy.copy(info['stdout_path'])
        stderr_path = copy.copy(info['stderr_path'])
        hints = copy.deepcopy(info['hints'])

        program_name = copy.deepcopy(args[0])
        if program_name.__class__ == list:
            new_args = args[0]
            program_name = new_args[-1]
            new_args.extend(args[1:])
            args = new_args

        self.check_subprocess_command(args)
        # launch the process and always pipe stdout and stderr because we
        # want to watch both streams, regardless of whether stdout should
        # be passed on to another process
        proc = subprocess.Popen(
            args,
            stdin = use_stdin,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE,
            preexec_fn = restore_sigpipe_handler,
            close_fds = True
        )
        pid = proc.pid
        self.popen_procs[pid] = proc

        self.running_procs.add(pid)
        self.proc_order.append(pid)
        self.proc_details[pid] = {
            'name': os.path.basename(program_name),
            'start_time': datetime.datetime.now(),
            'args': args,
            'pid': pid,
            'hints': hints
        }
        self.log("Launched %s as PID %d." % (' '.join(args), pid))
        sys.stderr.write("Launched %s as PID %d.\n" % (' '.join(args), pid))
        sys.stderr.flush()

        pipe = None
        if keep_stdout_open:
            pipe = os.pipe()

        self.copy_processes_for_pid[pid] = list()

        for which in ['stdout', 'stderr']:
            report_path = self.get_run().add_temporary_file("%s-report" % which)
            sink_path = stdout_path if which == 'stdout' else stderr_path
            listener_pid = self._do_launch_copy_process(
                proc.stdout if which == 'stdout' else proc.stderr,
                sink_path,
                report_path, pid, which,
                pipe if which == 'stdout' else None)

            self.copy_processes_for_pid[pid].append(listener_pid)
            self.copy_process_reports[listener_pid] = report_path

            if sink_path is not None:
                self.proc_details[listener_pid]['sink'] = os.path.basename(sink_path)
                self.proc_details[listener_pid]['sink_full_path'] = sink_path

        if keep_stdout_open:
            os.close(pipe[1])
            return pipe[0], pid
        else:
            return None, pid

    def _do_launch_copy_process(self, fin, fout_path, report_path, parent_pid, which, pipe):
        pid = os.fork()
        if pid == 0:

            def write_report_and_exit():
                try:
                    # write report
                    with open(report_path, 'w') as freport:
                        report = dict()
                        report['sha1'] = checksum.hexdigest()
                        report['tail'] = tail
                        report['length'] = length
                        report['lines'] = newline_count
                        freport.write(yaml.dump(report))
                except (IOError, LookupError) as e:
                    sys.stderr.write("Exception (%s): %s" % (type(e).__name__, sys.exc_info() ))
                    sys.stderr.write(traceback.format_exc())
                    sys.stderr.flush()
                    pass

                os._exit(0)

            def sigpipe_handler(signum, frame):
                write_report_and_exit()

            signal.signal(signal.SIGTERM, sigpipe_handler)
            signal.signal(signal.SIGINT, sigpipe_handler)
            signal.signal(signal.SIGPIPE, sigpipe_handler)
            os.setsid()
            if pipe is not None:
                os.close(pipe[0])
            fdout = None
            if fout_path is not None:
                fdout = os.open(fout_path, os.O_WRONLY|os.O_CREAT|os.O_TRUNC)

            checksum = hashlib.sha1()
            tail = ''
            length = 0
            newline_count = 0

            while True:
                block = fin.read(ProcessPool.COPY_BLOCK_SIZE)
                if len(block) == 0:
                    # fin reports EOF, let's call it a day
                    break

                # update checksum
                checksum.update(block)

                # update tail
                if len(block) >= ProcessPool.TAIL_LENGTH:
                    tail = block[-ProcessPool.TAIL_LENGTH:]
                else:
                    keep_length = ProcessPool.TAIL_LENGTH - len(block)
                    tail = tail[0:keep_length] + block

                # update length
                length += len(block)

                # update newline_count
                newline_count += block.count('\n')

                # write block to output file
                if fdout is not None:
                    bytes_written = os.write(fdout, block)
                    if bytes_written != len(block):
                        os._exit(1)

                # write block to pipe
                if pipe is not None:
                    bytes_written = os.write(pipe[1], block)
                    if bytes_written != len(block):
                        os._exit(2)

            # we're finished, close everything
            fin.close()
            if fdout is not None:
                os.close(fdout)
            if pipe is not None:
                os.close(pipe[1])

            write_report_and_exit()

            os._exit(0)
        else:
            self.running_procs.add(pid)
            self.proc_order.append(pid)
            self.proc_details[pid] = {
                'name': '[stream listener for %s of PID %d]' % (which, parent_pid),
                'start_time': datetime.datetime.now(),
                'pid': pid,
                'listener_for': [parent_pid, which]
            }
            self.log("Launched a copy process with PID %d to capture %s of PID %d." % (pid, which, parent_pid))
            if fout_path is not None:
                self.log("...which gets also redirected to %s" % fout_path)
            return pid

    def _wait(self):
        '''
        Wait for all processes to exit.
        '''
        self.log("Now launching process watcher and waiting for all child "
                 "processes to exit.")
        watcher_report_path = self.get_run()\
                                  .add_temporary_file('watcher-report')
        watcher_pid = self._launch_process_watcher(watcher_report_path)
        ProcessPool.process_watcher_pid = watcher_pid
        something_went_wrong = False
        pid = None
        while True:
            if len(self.running_procs) == 0:
                break
            try:
                # wait for the next child process to exit
                pid, exit_code_with_signal = os.wait()
                sys.stderr.write("PID: %s, Exit code: %s\n" %
                                 (pid, exit_code_with_signal))
                sys.stderr.flush()
                if pid == watcher_pid:
                    ProcessPool.process_watcher_pid = None
                    try:
                        with open(watcher_report_path, 'r') as f:
                            self.process_watcher_report = yaml.load(f)
                    except IOError as e:
                        sys.stderr.write("%s (%s): %s" %
                              (type(e).__name__, e.errno, e.strerror))
                        sys.stderr.write("Warning: Couldn't load watcher report from %s." %
                              watcher_report_path)
                        sys.stderr.flush()
                        pass
                    # the process watcher has terminated, which is cool, I guess
                    # (if it's the last child process, anyway)
                    continue

                try:
                    # remove pid from self.running_procs
                    self.running_procs.remove(pid)
                except KeyError as e:
                    sys.stderr.write("Key error(%s): %s" % (e.args, e.message))
                    sys.stderr.flush()
                    if pid != os.getpid():
                        sys.stderr.write("Note: Caught a process which we "
                                         "didn't know: %d.\n" % pid)
                        sys.stderr.flush()
                if pid in self.proc_details:
                    self.proc_details[pid]['end_time'] = datetime.datetime.now()

                signal_number = exit_code_with_signal & 255
                exit_code = exit_code_with_signal >> 8
                what_happened = "has exited with exit code %d" % exit_code
                if signal_number > 0:
                    what_happened = "has received signal %d" % signal_number
                    if signal_number in ProcessPool.SIGNAL_NAMES:
                        what_happened = ("has received %s (signal number %d)" %
                        (ProcessPool.SIGNAL_NAMES[signal_number], signal_number))

                if pid in self.proc_details:
                    self.log("%s (PID %d) %s." % (self.proc_details[pid]['name'],
                                                  pid, what_happened))
                else:
                    self.log("PID %d %s." % (pid, what_happened))

                if pid in self.proc_details:
                    if signal_number == 0:
                        self.proc_details[pid]['exit_code'] = exit_code
                    else:
                        self.proc_details[pid]['signal'] = signal_number
                        if signal_number in ProcessPool.SIGNAL_NAMES:
                            self.proc_details[pid]['signal_name'] = ProcessPool.SIGNAL_NAMES[signal_number]

                    # now kill it's preceding process, if this is from a pipeline
                    if 'use_stdin_of' in self.proc_details[pid]:
                        pidlist = list()
                        pidlist.append(self.copy_processes_for_pid\
                                       [self.proc_details[pid]['use_stdin_of']]\
                                       [0])
                        pidlist.append(self.copy_processes_for_pid\
                                       [self.proc_details[pid]['use_stdin_of']]\
                                       [1])
                        pidlist.append(self.proc_details[pid]['use_stdin_of'])
                        for kpid in pidlist:
                            self.log("Now killing %d, the predecessor of %d." %
                                     (kpid, pid))
                            self.ok_to_fail.add(kpid)
                            try:
                                os.kill(kpid, signal.SIGPIPE)
                            except OSError as e:
                                if e.errno == errno.ESRCH:
                                    self.log("Couldn't kill %d: no such "
                                             "process." % kpid)
                                    pass
                                else:
                                    raise

                    if pid in self.copy_process_reports:
                        report_path = self.copy_process_reports[pid]
                        if os.path.exists(report_path):
                            with open(report_path, 'r') as f:
                                report = yaml.load(f)

                            if report is not None:
                                self.proc_details[pid].update(report)

            except TimeoutException as e:
                sys.stderr.write("TimeoutException (%s): %s" % (e.args, e.message))
                sys.stderr.flush()
                self.log("Timeout, killing all child processes now.")
                ProcessPool.kill_all_child_processes()
            except OSError as e:
                if e.errno == errno.ECHILD:
                    # no more children running, we are done
                    sys.stderr.write("ProcessPool: There are no child "
                                     "processes left, exiting.\n")
                    sys.stderr.flush()
                    signal.alarm(0)
                    self.log("Cancelling timeout (if there was one), all "
                             "child processes have exited.")
                    break
                elif e.errno == errno.EINTR:
                    # a system call was interrupted, pfft.
                    pass
                else:
                    raise
            else:
                if exit_code_with_signal != 0:
                    if not pid in self.ok_to_fail:
                        # Oops, something went wrong. See what happens and
                        # terminate all child processes in a few seconds.
                        something_went_wrong = True
                        signal.signal(signal.SIGALRM, timeout_handler)
                        self.log("Terminating all children in %d seconds..." %
                                 ProcessPool.SIGTERM_TIMEOUT)
                        signal.alarm(ProcessPool.SIGTERM_TIMEOUT)

        # now wait for the watcher process, if it still exists
        try:
            os.waitpid(watcher_pid, 0)
            try:
                with open(watcher_report_path, 'r') as f:
                    self.process_watcher_report = yaml.load(f)
            except IOError as e:
                sys.stderr.write("%s (%s): %s" % (type(e).__name__, e.errno, e.strerror))
                sys.stderr.write("Warning: Couldn't load watcher report from %s." %
                      watcher_report_path)
                sys.stderr.flush()
                pass
        except OSError as e:
            if e.errno == errno.ESRCH:
                pass
            elif e.errno == errno.ECHILD:
                pass
            else:
                raise

        if something_went_wrong:
            log = "Pipeline crashed. PID: %s" % pid
            self.log(log)
            logger.error(log)
            sys.exit(1)

    def _launch_process_watcher(self, watcher_report_path):
        '''
        Launch the process watcher via fork. The process watcher repeatedly
        determines all child processes of the main process and determines their
        current and maximum CPU and RAM usage.
        Initially, this is done in short intervals (0.1 seconds), so that very
        short-lived processes can be watched but the frequency drops quickly so
        that after a while, child processes are only examined every 10 seconds.
        '''
        super_pid = os.getpid()

        watcher_pid = os.fork()
        if watcher_pid == 0:
            try:
                signal.signal(signal.SIGTERM, signal.SIG_DFL)
                signal.signal(signal.SIGINT, signal.SIG_IGN)
                called_cpu_stat_for_childpid = set()
                procs = {}
                procs[super_pid] = psutil.Process(super_pid)
                procs[os.getpid()] = psutil.Process(os.getpid())
                for pid in self.running_procs:
                    try:
                        procs[pid] = psutil.Process(pid)
                    except psutil.NoSuchProcess:
                        pass

                pid_list = copy.deepcopy(procs.keys())
                for pid in pid_list:
                    proc = procs[pid]
                    try:
                        cpu_percent = proc.cpu_percent(interval = None)
                        # add values for all children
                        if pid != super_pid:
                            for p in proc.children(recursive = True):
                                try:
                                    cpu_percent = p.cpu_percent(interval = None)
                                    called_cpu_stat_for_childpid.add(p.pid)
                                except psutil.NoSuchProcess:
                                    pass
                    except psutil.NoSuchProcess:
                        del procs[pid]

                time.sleep(0.1)

                iterations = 0
                delay = 0.5
                max_data = dict()
                while True:
                    pid_list = copy.deepcopy(procs.keys())
                    sum_data = dict()
                    for pid in pid_list:
                        proc = procs[pid]
                        try:
                            data = dict()
                            data['cpu_percent'] = proc.cpu_percent(interval = None)
                            data['memory_percent'] = proc.memory_percent()
                            memory_info = proc.memory_info()
                            data['rss'] = memory_info.rss
                            data['vms'] = memory_info.vms

                            # add values for all children
                            if pid != super_pid:
                                for p in proc.children(recursive = True):
                                    try:
                                        v = p.cpu_percent(interval = None)
                                        if p.pid in called_cpu_stat_for_childpid:
                                            data['cpu_percent'] += v
                                        called_cpu_stat_for_childpid.add(p.pid)
                                        data['memory_percent'] += p.memory_percent()
                                        memory_info = p.memory_info()
                                        data['rss'] += memory_info.rss
                                        data['vms'] += memory_info.vms
                                    except psutil.NoSuchProcess:
                                        pass

                            if not pid in max_data:
                                max_data[pid] = copy.deepcopy(data)
                            for k, v in data.items():
                                max_data[pid][k] = max(max_data[pid][k], v)

                            if len(sum_data) == 0:
                                for k, v in data.items():
                                    sum_data[k] = 0.0

                            for k, v in data.items():
                                sum_data[k] += v

                        except psutil.NoSuchProcess:
                            del procs[pid]

                    if not 'sum' in max_data:
                        max_data['sum'] = copy.deepcopy(sum_data)
                    for k, v in sum_data.items():
                        max_data['sum'][k] = max(max_data['sum'][k], v)

                    if len(procs) <= 2:
                        # there's nothing more to watch, write report and exit
                        # (now there's only the controlling python process and
                        # the process watcher itself
                        with open(watcher_report_path, 'w') as f:
                            report = dict()
                            report['max'] = max_data
                            f.write(yaml.dump(report, default_flow_style = False))

                        os._exit(0)

                    iterations += 1
                    if iterations == 10:
                        delay = 1
                    if iterations == 20:
                        delay = 2
                    if iterations == 30:
                        delay = 10
                    time.sleep(delay)
            except:
                sys.stderr.write("PID (%s) Process Watcher Exception: %s" %
                      (os.get_pid(), sys.exc_info()))
                sys.stderr.flush()
            finally:
                os._exit(0)
        else:
            return watcher_pid

    @classmethod
    def kill(cls):
        '''
        Kills all user-launched processes. After that, the remaining process
        will end and a report will be written.
        '''
        ProcessPool.process_pool_is_dead = True

        signal.signal(signal.SIGALRM, timeout_handler)
        signal.alarm(ProcessPool.SIGTERM_TIMEOUT)
        if ProcessPool.current_instance is not None:
            for pid in ProcessPool.current_instance.copy_processes_for_pid.keys():
                try:
                    os.kill(pid, signal.SIGTERM)
                except Exception as e:
                    sys.stderr.write("PID (%s) threw %s: %s" %
                          (pid, type(e).__name__, sys.exc_info()))
                    sys.stderr.write(traceback.format_exc())
                    sys.stderr.flush()
                    pass

    @classmethod
    def kill_all_child_processes(cls):
        '''
        Kill all child processes of this process by sending a SIGTERM to each of them.
        This includes all children which were not launched by this module, and
        their children etc.
        '''
        proc = psutil.Process(os.getpid())
        for p in proc.children(recursive = True):
            try:
                p.terminate()
            except psutil.NoSuchProcess:
                pass

