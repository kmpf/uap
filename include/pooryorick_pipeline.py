#!/usr/bin/python

'''
Copyright 2007, Nathan Coulter, pipeline.py@pooryorick.com

Licensed under the same Same terms as the PYTHON SOFTWARE FOUNDATION LICENSE VERSION 2, taking 'PSF' to mean the copyright holder of this work, and 'Python' to mean this work. See http://www.python.org/download/releases/2.5.1/license/
'''

import os
import signal
import sys

class PipelineError(Exception):
    pass

def pipeline(*cmds):
    '''
    Return the file descriptor of the pipe output and the pids of the pipe
    commands.

    # Like exec* functions, the first argument to a command contains the
    # name of the command itself.
    EDIT: This has been patched because I don't see why we should repeat
    ourselves at the beginning of each command. Also, commands are passed
    as lists now rather than tuples.

    The first command reads its input from sys.stdin, or from the file
    descriptor stored in kw['read'].  All error streams go to sys.stderr.

    # >>> cmd1 = ('echo', 'echo', 'batman')
    # >>> cmd2 = ('tr', 'tr', 'a', 'o')
    # >>> cmd3 = ('rev', 'rev')
    EDIT:
    >>> cmd1 = ['echo', 'batman']
    >>> cmd2 = ['tr', 'a', 'o']
    >>> cmd3 = ['rev']

    >>> fd, pids = pipeline(cmd1, cmd2, cmd3)
    >>> print os.read(fd, 1024),
    nomtob
   >>> print close(pids)
   [0, 0, 0]

    '''
    pids = []
    read = sys.stdin.fileno()
    for cmd in cmds:
        cmd.insert(1, cmd[0])
        read1, write1 = os.pipe()
        pid = os.fork()
        if not pid:
            os.dup2(read, 0)
            os.dup2(write1, 1)
            os.close(write1)
            try:
                os.execvp(cmd[0], cmd[1:])
            except:
                cls, inst = sys.exc_info()[:2]
                inst.args = (inst.args[0],
                    '%s:  %s' % (inst.args[1], ' '.join(cmd)))
                raise(cls, inst), inst.args
        pids.append(pid)
        os.close(write1)
        read = read1
    return read1, pids

def open(*cmds):
    '''
    like pipeline, but returns an open file instead of a file descriptor
    Example:

        >>> cmds = []
        >>> cmd = ('echo', 'echo', 'easy as pie')
        >>> cmds.append(cmd)
        >>> cmd = ('awk', 'awk', '{ print $1, $2, "rhubarb", $3 }')
        >>> cmds.append(cmd)
        >>> output, pids = open(*cmds)
        >>> print output.read(),
        easy as rhubarb pie


        >>> cmd1 = ('echo', 'echo', 'hello')
        >>> cmd2 = ('sed', 'sed', 's/hello/humbug/')
        >>> output, pids = open(cmd1, cmd2)
        >>> print output.read(),
        humbug

        >>> cmd1 = ('yes', 'yes', 'holy bat-pipe, Batman!')
        >>> cmd2 = ('sed', 'sed', 's/Batman/Smithers/')
        >>> output, pids = open(cmd1, cmd2)
        >>> print output.read(4),
        holy
        >>> print output.read(4)
         bat
        >>> print output.read(4)
        -pip

        >>> cmd4 = ('rev', 'rev')
        >>> from StringIO import StringIO
        >>> file_ = StringIO("hello")
        >>> output, pids = open(file_, cmd4)
        >>> print output.read(),
        olleh

    Also see 'example' functions in source code

    '''
    #stdout, pids = pipeline(*cmds)
    #return os.fdopen(stdout), pids

    if getattr(cmds[0], 'read', None):
        read, write = os.pipe()
        sys.stdin = os.fdopen(read)
        pid = os.fork()
        if not pid:
            sys.stdin.close()
            while 1:
                data = cmds[0].read(512)
                if not data: break
                os.write(write, data)
                os._exit(0)
        else:
            os.close(write)
            stdout, pids = pipeline(*cmds[1:])
    else:
        stdout, pids = pipeline(*cmds)
    return os.fdopen(stdout), pids


def read(*cmds):
    '''
    Returns the output a pipeline.  Raises an PipelineError if any component of
    the pipeline fails.

    Example:

    >>> cmd1 = ('echo', 'echo', 'green beans')
    >>> cmd2 = ('awk' , 'awk', '/beans/{ print $1 }')
    >>> cmd3 = ('rev' , 'rev')
    >>> cmd4 = ('bash', 'bash', '-c', 'read; echo $REPLY $REPLY')
    >>> print read(cmd1, cmd2, cmd3, cmd4),
    neerg neerg
    '''
    file_, pids = open(*cmds)
    out = file_.read()
    s = close(pids)
    if max(s) > 0:
        raise PipelineError, "return status:  %s" % s
    return out

def close(pids):
    '''
   Kills any running pipeline processes

   Returns the exit status of each process in the pipeline

    >>> cmd1 = ('echo', 'echo', 'hello')
    >>> cmd2 = ('tr', 'tr', 'e', 'a')
    >>> output, pids = open(cmd1, cmd2)
    >>> print output.read(),
    hallo
    >>> print close(pids)
    [0, 0]
    '''
    for pid in pids:
        os.kill(pid, signal.SIGTERM)
        os.kill(pid, signal.SIGKILL)
    return map(lambda x: os.WEXITSTATUS(os.waitpid(x,0)[1]), pids)


if __name__ == '__main__':

    def test1():
        line, pids = open(
            ('yes', 'yes', 'easy as py'),
            ('cat', 'cat'),
            ('tr', 'tr', 'a', 'A'),
            ('tr', 'tr', 'g', 'G'),
            ('awk', 'awk', '{ print $0 }'),
            ('sed', 'sed', 's/pipe/PIPE/g'),
            ('sed', 'sed', 's/import/IMPORT/g'),
        )

        print line.read(3)
        print line.read(3)
        line.close()

        #for i in range(100):
        #    print line.read(3),


    def test2():
        cmd2 = ('tr', 'tr', 'a', 'o')
        cmd3 = ('rev', 'rev')
        read, write = os.pipe()
        fd, pids = pipeline(cmd2, cmd3)
        pid = os.fork()
        if not pid:
            os.close(read)
            os.write(write, 'tack')
        else:
            os.close(write)
            print os.read(fd, 1024)


    import doctest
    doctest.testmod()
    #test1()
    #test2()
