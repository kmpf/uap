import sys
import yaml

from ..abstract_step import *
from .. import process_pool

logger = logging.getLogger('uap_logger')

# Taken from:
# http://stackoverflow.com/questions/2023608/check-what-files-are-open-in-python
#import __builtin__
#openfiles = set()
#oldfile = __builtin__.file
#class newfile(oldfile):
#    def __init__(self, *args):
#        self.x = args[0]
#        print "### OPENING %s ###" % str(self.x)
#        oldfile.__init__(self, *args)
#        openfiles.add(self)
#        printOpenFiles()
#
#    def close(self):
#        print "### CLOSING %s ###" % str(self.x)
#        oldfile.close(self)
#        openfiles.remove(self)
#oldopen = __builtin__.open
#def newopen(*args):
#    return newfile(*args)
#__builtin__.file = newfile
#__builtin__.open = newopen
#
#def printOpenFiles():
#    print "### %d OPEN FILES: [%s]" % (len(openfiles), ", ".join(f.x for f in openfiles))
#    print "### %d OPEN FILES" % len(openfiles)

class FastxQualityStats(AbstractStep):

    def __init__(self, pipeline):
        super(FastxQualityStats, self).__init__(pipeline)
        
        self.set_cores(4)
        
        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/first_read_quality_stats')
        self.add_connection('out/second_read_quality_stats')

        self.add_option('new_output_format', bool, default= True, optional=True)
        self.add_option('quality', int, default=33, optional=True)
        
        self.require_tool('cat')
        self.require_tool('dd')
        self.require_tool('mkfifo')
        self.require_tool('fastx_quality_stats')
        self.require_tool('pigz')

    def runs(self, run_ids_connections_files):
        # get a list of all read files we have to count
        sample_input_paths_dict = dict()
        reads_counts_files = dict()
        read_files = list()

        read_types = {'first_read': '_R1', 'second_read': '_R2'}
        for run_id in list(run_ids_connections_files.keys()):
            with self.declare_run(run_id) as run:
                for read in read_types:
                    connection = 'in/%s' % read
                    input_paths = run_ids_connections_files[run_id][connection]
                    
                    # Check for empty connections
                    if input_paths == [None]:
                        run.add_empty_output_connection(
                            "%s_quality_stats" % read)
                    else:
                        temp_fifos = list()
                        exec_group = run.new_exec_group()
                        for input_path in input_paths:
                            temp_fifo = run.add_temporary_file(
                                "fifo-%s" %
                                os.path.basename(input_path) )
                            temp_fifos.append(temp_fifo)
                            mkfifo = [self.get_tool('mkfifo'), temp_fifo]
                            exec_group.add_command(mkfifo)

                            # 2. Output files to fifo
                            if input_path.endswith('fastq.gz'):
                                with exec_group.add_pipeline() as unzip_pipe:
                                    # 2.1 command: Read file in 4MB chunks
                                    dd_in = [self.get_tool('dd'),
                                           'ibs=4M',
                                           'if=%s' % input_path]
                                    # 2.2 command: Uncompress file to fifo
                                    pigz = [self.get_tool('pigz'),
                                            '--decompress',
                                            '--stdout']
                                    # 2.3 Write file in 4MB chunks to fifo
                                    dd_out = [self.get_tool('dd'),
                                              'obs=4M',
                                              'of=%s' % temp_fifo]
                                
                                    unzip_pipe.add_command(dd_in)
                                    unzip_pipe.add_command(pigz)
                                    unzip_pipe.add_command(dd_out)
                            elif input_path.endswith('fastq'):
                                # 2.1 command: Read file in 4MB chunks and
                                #              write to fifo in 4MB chunks
                                dd_in = [self.get_tool('dd'),
                                         'bs=4M',
                                         'if=%s' % input_path,
                                         'of=%s' % temp_fifo]
                                exec_group.add_command(dd_in)
                            else:
                                raise Exception("File %s does not end with "
                                                    "any expected suffix ("
                                                    "fastq.gz or fastq). Please "
                                                    "fix that issue.")

                        # 3. Read data from fifos and check quality stats
                        with exec_group.add_pipeline() as fastx_pipe:
                            # 3.1 command: Read from ALL fifos
                            cat = [self.get_tool('cat')]
                            cat.extend(temp_fifos)
                            # 3.2 command: Compute quality statistics
                            fastx_qs_file = run.add_output_file(
                                "%s_quality_stats" % read,
                                "%s%s.fastq.quality.tsv" %
                                (run_id, read_types[read]),
                                input_paths)
                            fastx_qs = [self.get_tool('fastx_quality_stats'),
                                        '-Q', '%s' % self.get_option('quality')]
                            if self.get_option('new_output_format'):
                                fastx_qs.append('-N')

                            fastx_pipe.add_command(cat)
                            fastx_pipe.add_command(fastx_qs,
                                                   stdout_path = fastx_qs_file)
