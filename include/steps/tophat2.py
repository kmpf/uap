import sys
from abstract_step import *
import glob
import misc
import process_pool
import yaml


class TopHat2(AbstractStep):
    '''
    TopHat is a fast splice junction mapper for RNA-Seq reads.
    It aligns RNA-Seq reads to mammalian-sized
    genomes using the ultra high-throughput short read aligner Bowtie
    , and then analyzes the mapping results to identify splice junctions between exons.

    http://tophat.cbcb.umd.edu/
    
    typical command line::

        tophat [options]* <index_base> <reads1_1[,...,readsN_1]> [reads1_2,...readsN_2]


    '''
    
    def __init__(self, pipeline):
        super(TopHat2, self).__init__(pipeline)
        self.set_cores(6)
        
        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/alignments')
        self.add_connection('out/unmapped')
        self.add_connection('out/insertions')
        self.add_connection('out/deletions')
        self.add_connection('out/junctions')
        self.add_connection('out/misc_logs')
        self.add_connection('out/log_stderr')
        self.add_connection('out/prep_reads')
        self.add_connection('out/align_summary')

        self.require_tool('cat4m')
        self.require_tool('pigz')
        self.require_tool('bowtie2')
        self.require_tool('tophat2')

#        self.add_option('genome', str)
        self.add_option('index', str)

        self.add_option('library_type', str, choices = [
            'fr-unstranded', 'fr-firststrand', 'fr-secondstrand'])
        
        # TODO: remove swap_reads option
        self.add_option('swap_reads', bool, default = False)

    def declare_runs(self):
        ### make sure files are available
        if not os.path.exists(self.get_option('index') + '.1.bt2'):
            raise StandardError("Could not find index file: %s.*" % self.get_option('index'))

        found_files = dict()
        read_types = {'first_read': '-R1', 'second_read': '-R2'}
        paired_end_info = dict()

        for read in read_types.keys():
            for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/%s' % read):
                if input_paths != [None]:
                    paired_end_info[run_id] = self.find_upstream_info_for_input_paths(input_paths, 'paired_end')

                    # save information per file in found_files
                    if not run_id in found_files:
                        found_files[run_id] = dict()

                    if not read in found_files[run_id]:
                        found_files[run_id][read] = list()
                    # Check if we get exactly one input file
                    if len(input_paths) != 1:
                        raise StandardError("Expected one input file.")
                    found_files[run_id][read].extend(input_paths)


        for run_id in found_files.keys():
            with self.declare_run(run_id) as run:
                run.new_exec_group()
                run.add_private_info('paired_end', paired_end_info[run_id])
                input_paths = list()
                for read in found_files[run_id].keys():
                    run.add_private_info(read, found_files[run_id][read])
                    input_paths.extend(found_files[run_id][read])

                run.add_output_file(
                    'alignments', 
                    '%s-tophat2-accepted.bam' % run_id,input_paths)
                run.add_output_file(
                    'unmapped', 
                    '%s-tophat2-unmapped.bam' % run_id, input_paths)
                run.add_output_file(
                    'log_stderr', 
                    '%s-tophat2-log_stderr.txt' % run_id, input_paths)
                run.add_output_file(
                    'misc_logs', 
                    '%s-tophat2-misc_logs.yaml' % run_id, input_paths)
                run.add_output_file(
                    'insertions', 
                    '%s-tophat2-insertions.bed' % run_id, input_paths)
                run.add_output_file(
                    'junctions', 
                    '%s-tophat2-junctions.bed' % run_id, input_paths)
                run.add_output_file(
                    'deletions', 
                    '%s-tophat2-deletions.bed' % run_id, input_paths)
                run.add_output_file(
                    'prep_reads', 
                    '%s-tophat2-prep_reads.info' % run_id, input_paths)
                run.add_output_file(
                    'align_summary', 
                    '%s-tophat2-align_summary.txt' % run_id, input_paths)

    def execute(self, run_id, run):
        read_types = {'first_read': '-R1', 'second_read': '-R2'}
        is_paired_end = run.get_private_info('paired_end')
        first_read_path = run.get_private_info('first_read')[0]
        second_read_path = None
        if is_paired_end:
            second_read_path = run.get_private_info('second_read')[0]

        tophat2_out_path = self.get_temporary_path('tophat-out')

        with process_pool.ProcessPool(self) as pool:
            with pool.Pipeline(pool) as pipeline:

                cores = str(self.get_cores())
                tophat2 = [
                    self.get_tool('tophat2'),
                    '--library-type', self.get_option('library_type'),
                    '--output-dir', tophat2_out_path,
                    '-p', cores, self.get_option('index'), 
                    first_read_path
                    ]
                
                if is_paired_end:
                    tophat2.append(second_read_path)

                log_stderr = run.get_single_output_file_for_annotation(
                    'log_stderr')
                pipeline.append(tophat2, stderr_path = log_stderr)
                
        tophat2_generic_names = [ 'accepted_hits.bam', 'unmapped.bam', 'insertions.bed', 'deletions.bed', 'junctions.bed', 'prep_reads.info', 'align_summary.txt']
        for i  in  tophat2_generic_names:
            target_name = i.split('.')[0]
            if (target_name == 'accepted_hits'):
                target_name = 'alignments'
            try:
                os.rename(os.path.join(tophat2_out_path, i), run.get_single_output_file_for_annotation(target_name))
            except OSError:
                raise StandardError( '\n\t Something went wrong while renaming tophat2 output files: \
                                      \n\t old_dir       : {0} \
                                      \n\t old_basename  : {1} \
                                      \n\t key_for_target: {2} \
                                      \n\t target_path   : {3} '.format(tophat2_out_path, i, target_name, run.get_single_output_file_for_annotation(target_name)))

        misc_logs = run.get_single_output_file_for_annotation('misc_logs')
        with open(misc_logs, 'w') as f:
            logs = {}
            for path in glob.glob(os.path.join(tophat2_out_path, 'logs', '*')):
                logs[os.path.basename(path)] = open(path).read()
                os.unlink(path)
            try:
                os.rmdir(os.path.join(tophat2_out_path, 'logs'))
            except OSError:
                pass
            f.write(yaml.dump(logs, default_flow_style = False))

        try:
            os.rmdir(tophat2_out_path)
        except OSError:
            pass
