import sys
from abstract_step import *
import process_pool
import re
import yaml


class BedToBedgraph(AbstractStep):

    def __init__(self, pipeline):
        super(BedToBedgraph, self).__init__(pipeline)
        
        self.set_cores(4)
        
        self.add_connection('in/alignments')
        self.add_connection('out/tracks')
        
        self.require_tool('cat4m')
        self.require_tool('bedtools')

        self.add_option('genome', str, optional=False)
        self.add_option('strand-specific', bool, optional=False)

    def declare_runs(self):
        is_strand_specific = self.get_option('strand-specific')
        for run_id, input_paths in self.get_run_ids_and_input_files_for_connection('in/alignments'):
            with self.declare_run(run_id) as run:
                
                if len(input_paths) != 1:
                    raise StandardError("Expected exactly one alignment file.")
                if  input_paths[0][-4:] != '.bed':
                    raise StandardError("%s file suffix is not '.bed'. " +
                                        "Please provide a BED file" % input_paths[0])

                run.add_private_info('in-bed', input_paths[0])
                if is_strand_specific:
                    plus_file = os.path.basename(input_paths[0])[:-4] + '.plus.bedgraph'
                    minus_file = os.path.basename(input_paths[0])[:-4] + '.minus.bedgraph'
                    run.add_output_file('tracks', plus_file, input_paths)
                    run.add_output_file('tracks', minus_file, input_paths)
                    bg_files = [plus_file, minus_file]
                    run.add_private_info('bg-files', bg_files)

                    plus_trackopts = '\"name=\"%s\"+ description=\"%s plus strand\"\"' % (run_id, run_id)
                    minus_trackopts = '\"name=\"%s\"- description=\"%s minus strand\"\"' % (run_id, run_id)
                    bg_trackopts = [plus_trackopts, minus_trackopts]
                    run.add_private_info('bg-trackopts', bg_trackopts)
                else:
                    bg_file = os.path.basename(input_paths[0])[:-4] + '.bedgraph'
                    run.add_output_file('tracks', bg_file, input_paths)
                    run.add_private_info('bg-files', [bg_file])

                    bg_trackopts = ['\"name=\"%s\" description=\"%s\"\"' % (run_id, run_id)]
                    run.add_private_info('bg-trackopts', bg_trackopts)

    def execute(self, run_id, run):
        is_strand_specific = self.get_option('strand-specific')
        bed_path = run.get_private_info('in-bed')
        bg_files = run.get_private_info('bg-files')
        bg_trackopts = run.get_private_info('bg-trackopts')
        track_info = dict()
        if is_strand_specific:
            if len(bg_files) == 2 and len(bg_trackopts) == 2:
                files = misc.assign_strings(bg_files, ['plus','minus'])
                trackopts = misc.assign_strings(bg_trackopts, ['plus','minus'])
                for strand in ['plus', 'minus']:
                    track_info[files[strand]] = trackopts[strand]
            else:
                raise StandardError("Expected two output files and track infos, " +
                                    "but instead got this output files %s and track " +
                                    "infos %s" % (bg_files, bg_trackopts))
        else:
            if len(bg_files) == 1 and len(bg_trackopts) == 1:
                track_info[bg_files[0]] = bg_trackopts[0]
            else:
                raise StandardError("Expected two output files and track infos, " +
                                    "but instead got this output files %s and track " +
                                    "infos %s" % (bg_files, bg_trackopts))

        for track_path, trackopts in track_info.items():
            with process_pool.ProcessPool(self) as pool:
                with pool.Pipeline(pool) as pipeline:
                    cat4m_in = [self.get_tool('cat4m'), bed_path]
                    bedtools = [self.get_tool('bedtools'), 'genomecov',
                                '-trackline',
                                '-trackopts', trackopts,
                                '-bg', '-split'
                                ]
                    strand = None
                    if is_strand_specific:
                        strand = misc.assign_string(track_path, ['plus','minus'])
                        if strand == 'plus':
                            bedtools.extend(['-strand', '+'])
                        elif strand == 'minus':
                            bedtools.extend(['-strand', '-'])
                        
                    bedtools.extend(['-g', self.get_option('genome'), '-i', 'stdin'])
#
                    cat4m_out = [self.get_tool('cat4m'), '-']
                    out_path = None
                    if is_strand_specific:
                        out_paths = run.get_output_files_for_annotation_and_tags('tracks', ['plus', 'minus'])
                        out_path = out_paths[strand]
                    else:
                        out_path = run.get_single_output_file_for_annotation('tracks')

                    pipeline.append(cat4m_in)
                    pipeline.append(bedtools, stdout_path = out_path)

