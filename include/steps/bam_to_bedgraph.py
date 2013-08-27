import sys
from abstract_step import *
import process_pool
import re
import yaml


class BamToBedgraph(AbstractStep):

    def __init__(self, pipeline):
        super(BamToBedgraph, self).__init__(pipeline)
        
        self.set_cores(4)
        
        self.add_connection('in/alignments')
        self.add_connection('out/tracks')
        
        self.require_tool('cat4m')
        self.require_tool('bedtools')

        self.add_option('strand-specific', bool, optional=False)

    def declare_runs(self):
        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/alignments'):
            with self.declare_run(run_id) as run:
                
                if len(input_paths) != 1:
                    raise StandardError("Expected exactly one alignment file.")
                if  input_paths[0][-4:] != '.bam':
                    raise StandardError("%s file suffix is not '.bam'. " +
                                        "Please provide a BAM file" % input_paths[0])

                if self.get_option('strand-specific'):
                    plus_file = os.path.basename(input_paths[0])[:-4] + '.plus.bedgraph'
                    minus_file = os.path.basename(input_paths[0])[:-4] + '.minus.bedgraph'
                    run.add_output_file('tracks', plus_file, input_paths)
                    run.add_output_file('tracks', minus_file, input_paths)
                    bg_files = [plus_file, minus_file]
                    run.add_private_info('bg-files', bg_files)

                    plus_trackopts = 'name=%s+ description=\"%s plus strand\"' % (run_id, run_id)
                    minus_trackopts = 'name=%s- description=\"%s minus strand\"' % (run_id, run_id)
                    bg_trackopts = [plus_trackopts, minus_trackopts]
                    run.add_private_info('bg-trackopts', bg_trackopts)
                else:
                    bg_file = os.path.basename(input_paths[0])[:-4] + '.bedgraph'
                    run.add_output_file('tracks', bg_file, input_paths)
                    run.add_private_info('bg-files', [bg_file])

                    bg_trackopts = 'name=%s description=\"%s\"' % (run_id, run_id)
                    run.add_private_info('bg-trackopts', bg_trackopts)


    def execute(self, run_id, run):
        bam_path = run.get_private_info('in-bam')
        bg_files = run.get_private_info('bg-files')
        bg_trackopts = run.get_private_info('bg_trackopts')
        track_info = dict()
        if len(bg_files) == 2 and len(bg_trackopts) == 2:
            files = misc.assign_strings(bg_files, ['plus','minus'])
            trackopts = misc.assign_strings(bg_trackopts, ['plus','minus'])
            for strand in ['plus', 'minus']:
                track_info[files[strand]] = trackopts[strand]
        elif len(bg_files) == 1 and len(bg_trackopts) == 1:
            track_info[bg_files[0]] = bg_trackopts
        else:
            raise StandardError("These track files %s and these trackopts " +
                                "%s do not match our criteria" % (bg_files, bg_trackopts))

        for track_path, trackopts in track_info.items():
            with process_pool.ProcessPool(self) as pool:
                with pool.Pipeline(pool) as pipeline:
                    cat4m_in = [self.get_tool('cat4m'), bam_path]
                    bedtools = [self.get_tool('bedtools'), 'genomecov',
                                '-trackline',
                                '-trackopts', trackopts,
                                '-bg', '-split'
                                ]
                    strand = misc.assign_string(track_path, ['plus','minus'])
                    if strand == 'plus':
                        bedtools.extend(['-strand', '+'])
                    elif strand == 'minus':
                        bedtools.extend(['-strand', '-'])
                    bedtools.extend(['-ibam', 'stdin'])

                    cat4m_out = [self.get_tool('cat4m'), '-']
                
                    pipeline.append(cat4m_in)
                    pipeline.append(bedtools)
                    pipeline.append(cat4m_out, stdout_path = track_path)
