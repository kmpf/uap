import sys
from abstract_step import *
import glob
import misc
import process_pool
import yaml
import os

from logging import getLogger

logger = getLogger('uap_logger')


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

        # all .gft assemblies from all samples that have been produced with
        # cufflinks
        self.add_connection('in/features')
        # merged assembly 'merged.gft'
        self.add_connection('out/features')  # merged.gtf
        self.add_connection('out/assemblies')  # input assemblies txt file
        self.add_connection('out/log_stderr')
        self.add_connection('out/run_log')

        self.require_tool('cuffmerge')
        self.require_tool('printf')
        self.require_tool('mkdir')
        self.require_tool('mv')

        # output dir (-o option) is fixed by out connection ?!
        self.add_option('run_id', str, optional=True,
                        description='An arbitrary name of the new '
                        'run (which is a merge of all samples).',
                        default='magic')
        self.add_option(
            'ref-gtf',
            str,
            optional=True,
            description='A "reference" annotation GTF. The input assemblies are merged together with the reference GTF and included in the final output.')
        self.add_option(
            'ref-sequence',
            str,
            optional=True,
            description='This argument should point to the genomic DNA sequences for the reference. If a directory, it should contain one fasta file per contig. If a multifasta file, all contigs should be present.')
        self.add_option(
            'num-threads',
            int,
            optional=True,
            description='Use this many threads to merge assemblies.',
            default=self.get_cores())

    def runs(self, run_ids_connections_files):

        # compile list of options
        options = ['ref-gtf', 'ref-sequence', 'num-threads']
        file_options = ['ref-gtf', 'ref-sequence']

        set_options = [option for option in options if
                       self.is_option_set_in_config(option)]

        option_list = list()
        for option in set_options:
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    option_list.append('--%s' % option)
            else:
                value = str(self.get_option(option))
                if option in file_options:
                    value = os.path.abspath(value)
                option_list.append('--%s' % option)
                option_list.append(value)

        # get all paths to the cufflinks assemblies from each sample
        cufflinks_sample_gtf = []
        for run_id in run_ids_connections_files.keys():
            cufflinks_sample_gtf.append(
                run_ids_connections_files[run_id]['in/features'][0])

        run_id = self.get_option('run_id')
        with self.declare_run(run_id) as run:

            # create the filename of the assemblies.txt file
            assemblies = [
                self.get_tool('printf'),
                '\n'.join(cufflinks_sample_gtf)]

            assemblies_file = run.add_output_file(
                'assemblies', '%s-cuffmerge-assemblies.txt' %
                run_id, cufflinks_sample_gtf)

            # 1. create assemblies file
            with run.new_exec_group() as as_exec_group:
                as_exec_group.add_command(
                    assemblies, stdout_path=assemblies_file)

            features_file = run.add_output_file(
                'features', '%s-cuffmerge-merged.gtf' %
                run_id, cufflinks_sample_gtf)
            run_log_file = run.add_output_file(
                'run_log', '%s-cuffmerge-run.log' %
                run_id, cufflinks_sample_gtf)
            log_err_file = run.add_output_file(
                'log_stderr', '%s-cuffmerge-log_stderr.txt' %
                run_id, cufflinks_sample_gtf)

            cuffmerge_out_path = run.add_temporary_directory('cuffmerge-out')
            cuffmerge = [self.get_tool('cuffmerge'), '-o', cuffmerge_out_path]
            cuffmerge.extend(option_list)

            # 2. Create temporary directory for cufflinks in- and output
            with run.new_exec_group() as dir_exec_group:
                mkdir = [self.get_tool('mkdir'), cuffmerge_out_path]
                dir_exec_group.add_command(mkdir)

            # 3. run cuffmerge
            with run.new_exec_group() as cm_exec_group:
                cuffmerge.append(assemblies_file)
                cm_exec_group.add_command(cuffmerge, stderr_path=log_err_file)

            result_files = {
                'merged.gtf': features_file,
                'logs/run.log': run_log_file,
            }

            # 4. mv output files from temp dir to final location
            with run.new_exec_group() as mv_exec_group:
                for orig, dest_path in result_files.items():
                    orig_path = os.path.join(cuffmerge_out_path, orig)
                    mv = [self.get_tool('mv'), orig_path, dest_path]
                    mv_exec_group.add_command(mv)

        # print assemblies to assemblies.txt (this is a run or so..)

#        sys.exit()
#
#        for run_id in run_ids_connections_files.keys():
#
#            with self.declare_run(run_id) as run:
#                input_paths = run_ids_connections_files[run_id]['in/features']
#                temp_dir = run.add_temporary_directory('cuffmerge-out')
#
#
#                cuffmerge = [self.get_tool('cuffmerge'), '-o', temp_dir]
#                cuffmerge.extend(option_list)
#
#
#                # 1. Create temporary directory for cufflinks in- and output
#                with run.new_exec_group() as exec_group:
#                    mkdir = [self.get_tool('mkdir'), temp_dir]
#                    exec_group.add_command(mkdir)
#
#                # 2. write assemblies file and append it to the cuffmerge cmd
#                with run.new_exec_group() as exec_group:
#                    assemblies = [self.get_tool('printf'), '\n'.join(input_paths)]
#                    assemblies_file = run.add_output_file('assemblies',
#                                                          '%s-cuffmerge-assemblies.txt' % run_id,
#                                                          input_paths)
#                    exec_group.add_command(assemblies,
#                                           stdout_path = assemblies_file)
#
#                    cuffmerge.append(assemblies_file)
#
#                result_files = {
#                    'merged.gtf': run.add_output_file('features',
#                                                      '%s-merged.gtf' % run_id,
#                                                      input_paths),
#                    'logs/run.log': run.add_output_file('run_log',
#                                                        '%s-run.log' % run_id,
#                                                        input_paths)
#                }
#
#                # 3. Execute cuffmerge
#                with run.new_exec_group() as exec_group:
#
# print(cuffmerge)
#
#                    exec_group.add_command(cuffmerge,
#                                           stderr_path = run.add_output_file('log_stderr',
#                                                                          '%s-cuffmerge-err.txt' % run_id,
#                                                                          input_paths
#                                                                          )
#                                           )
#
#                # 4. move files from temp-dir to usual uap-temp-dir, rename
#                with run.new_exec_group() as mv_exec_group:
#                    for orig, dest_path in result_files.items():
#                     # 3. Rename files
#                     orig_path = os.path.join(temp_dir, orig)
#                     mv = [self.get_tool('mv'), orig_path, dest_path]
#                     mv_exec_group.add_command(mv)
