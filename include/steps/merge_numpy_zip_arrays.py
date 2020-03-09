from uaperrors import StepError
import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger = getLogger('uap_logger')


class mergeNumpyZipArrays(AbstractStep):
    '''
    This step can be used to concatenate multiple zipped Numpy arrays which are
    the output of deepTools multiBamSummary subcommand.

    Usage example::

        merge_numpy_arrays.py [-h] [file-1.npz file-2.npz ... file-n.npz]  files

    '''

    def __init__(self, pipeline):
        super(mergeNumpyZipArrays, self).__init__(pipeline)

        self.set_cores(10)

        self.add_connection('in/read-coverage')
        self.add_connection('out/read-coverage')

        self.require_tool('merge_numpy_arrays.py')

    def runs(self, run_ids_connections_files):
        for run_id in run_ids_connections_files.keys():
            # Collect input_paths and labels for multiBamSummary
            input_paths = run_ids_connections_files[run_id]['in/alignments']
            labels = list()
            for f in input_paths:
                if not f.endswith(".bam"):
                    raise StepError(self, "Not a BAM file: %s" % f)
                if len(input_paths) > 1:
                    labels.append("%s-%s" % (run_id, input_paths.index(f)))
                else:
                    labels.append(run_id)

            with self.declare_run(run_id) as run:
                # Let's compile the command
                with run.new_exec_group() as multi_bam_summary_eg:
                    # 1. multiBamSummary command
                    # TODO: where is subcommand coming from???
                    multi_bam_summary = [
                        self.get_tool('multiBamSummary'), subcommand]
                    # Append list of input BAM files
                    multi_bam_summary.append('--bamfiles')
                    multi_bam_summary.extend(input_paths)
                    # Append name of the output file
                    multi_bam_summary.append('--outFileName')
                    multi_bam_summary.append(
                        run.add_output_file(
                            'read-coverage',
                            '%s.npz' % run_id,
                            input_paths
                        )
                    )
                    # Append list of BED files for BED-file subcommand
                    if subcommand == "BED-file":
                        multi_bam_summary.append('--BED')
                        multi_bam_summary.extend(self.get_option('bed-file'))
                    # Append list of labels
                    multi_bam_summary.append('--labels')
                    multi_bam_summary.extend(labels)
                    # Append number of processors
                    multi_bam_summary.extend(['--numberOfProcessors',
                                              str(self.get_cores())])
                    # Append list of options
                    # TODO: optionlist dont exists!!!
                    multi_bam_summary.extend(option_list)

                    # Add multiBamSummary to execution group
                    multi_bam_summary_eg.add_command(multi_bam_summary)
