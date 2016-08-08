import sys
from abstract_step import *
import glob
import misc
import process_pool
import yaml
import os

from logging import getLogger

logger=getLogger('uap_logger')

class CuffMerge(AbstractStep):

    '''
    CuffMerge is part of the 'Cufflinks suite of tools' for
    differential expr. analysis of RNA-Seq data and their
    visualisation. This step applies the cuffmerge tool which merges
    several Cufflinks assemblies. For details on cuffmerge we refer to
    the author's webpage:

    http://cole-trapnell-lab.github.io/cufflinks/cuffmerge/

    '''

    def __init__(self, pipeline):
        super(CuffMerge, self).__init__(pipeline)

        self.set_cores(6)

        self.add_connection('in/features') # ??? 2 assemplies per sample (TH + SM)
        self.add_connection('out/features') # merged.gt
        self.add_connection('out/assemblies_txt') # input assemblies txt file
        self.add_connection('out/log_stderr')
        self.add_connection('out/run_log')

        self.require_tool('cuffmerge')
        self.require_tool('printf')
        self.require_tool('mkdir')
        self.require_tool('mv')

        # output dir (-o option) is fixed by out connection ?!
        self.add_option('ref-gtf', str, optional=True,
                        description='A "reference" annotation GTF. The input assemblies are merged together with the reference GTF and included in the final output.')
        self.add_option('ref-sequence', str, optional=True,
                        description='This argument should point to the genomic DNA sequences for the reference. If a directory, it should contain one fasta file per contig. If a multifasta file, all contigs should be present.')
        self.add_option('num-threads', int, optional=True,
                        description='Use this many threads to merge assemblies.',
                        default = self.get_cores())

    def runs(self, run_ids_connections_files):
        
        # compile list of options
        options=['ref-gtf', 'ref-sequence', 'num-threads']

        set_options = [option for option in options if \
                       self.is_option_set_in_config(option)]

        option_list = list()
        for option in set_options:
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    option_list.append('--%s' % option)
            else:
                option_list.append('--%s' % option)
                option_list.append(str(self.get_option(option)))

        # number of threads is set to number of cores for this step.
        option_list.append('--num-threads %d' % int(self.get_cores()))

        for run_id in run_ids_connections_files.keys():

            with self.declare_run(run_id) as run:
                input_paths = run_ids_connections_files[run_id]['in/features']
                temp_dir = run.add_temporary_directory('cuffmerge-out')

                cuffmerge = [self.get_tool('cuffmerge'),
                             '-o', temp_dir]
                cuffmerge.extend(option_list)

                result_files = {
                    'merged.gtf': run.add_output_file(
                        'features',
                        '%s-merged.gtf' % run_id,
                        input_paths
                        ),
                    'logs/run.log': run.add_output_file(
                        'run_log',
                        '%s-run.log' % run_id,
                        input_paths
                        ),
                    'assemblies.txt': run.add_output_file(
                        'assemblies_txt',
                        '%s-cuffmerge-assemblies.txt' % run_id,
                        input_paths
                        )
                    }
                
                # 1. Create temporary directory for cufflinks in- and output
                with run.new_exec_group() as exec_group:
                    mkdir = [self.get_tool('mkdir'), temp_dir]
                    exec_group.add_command(mkdir)
                    
                # 2. Create Text file "manifest" with a list (one
                # per line) of GTF files that you'd like to merge
                # together into a single GTF file.
                with run.new_exec_group() as exec_group:

                    manifest = [self.get_tool('printf'), '\n'.join(input_paths)]
                    exec_group.add_command(manifest, stdout_path = run.add_output_file('assemblies.txt',
                                                                                       '%s-cuffmerge-assemblies.txt' % run_id,
                                                                                       input_paths)
                                           )
                    cufflinks.append(manifest)
                    
                # 3. Execute cuffmerge
                with run.new_exec_group() as exec_group:

                    exec_group.add_command(cuffmerge,
                                           stderr_path = run.add_output_file('log_stderr',
                                                                          '%s-cuffmerge-err.txt' % run_id,
                                                                          input_paths
                                                                          )
                                           )

                # 4. move files from temp-dir to usual uap-temp-dir, rename
                with run.new_exec_group() as exec_group:
                    for orig, dest_path in result_files.iteritems():
                     # 3. Rename files 
                     orig_path = os.path.join(temp_dir, orig)
                     mv = [self.get_tool('mv'), orig_path, dest_path]
                     mv_exec_group.add_command(mv)

                    
