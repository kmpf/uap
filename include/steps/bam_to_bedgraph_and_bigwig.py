from uaperrors import StepError
import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class BamToBedgraphAndBigwig(AbstractStep):

    def __init__(self, pipeline):
        super(BamToBedgraphAndBigwig, self).__init__(pipeline)

        self.set_cores(8)

        self.add_connection('in/alignments',
                            constraints = {'min_files_per_run': 1,
                                           'max_files_per_run': 1})
        self.add_connection('out/bedgraph')
        self.add_connection('out/bigwig')

        self.require_tool('bedtools')
        self.require_tool('sort')
        self.require_tool('bedGraphToBigWig')

        # General Options
        self.add_option('chromosome-sizes', str, optional = False)
        self.add_option('temp-sort-dir', str, optional = True)

    def runs(self, run_ids_connections_files):
        # Check if chromosome sizes points to a real file
        if not os.path.isfile(self.get_option('chromosome-sizes')):
            raise StepError(self, "Value for option 'chromosome-sizes' is not a "
                         "file: %s" % self.get_option('chromosome-sizes'))
        if self.get_option('temp-sort-dir') and \
           not os.path.isdir(self.get_option('temp-sort-dir')):
            raise StepError(self, "Value for option 'temp-sort-dir' is not a "
                         "directory: %s" % self.get_option('temp-sort-dir'))
        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                # Collect input paths
                input_paths = run_ids_connections_files[run_id]["in/alignments"]
                # Handle special condition e.g. no input files
                if input_paths == [None]:
                    run.add_empty_output_connection("alignments")
                # Complain if necessary
                elif len(input_paths) != 1:
                    raise StepError(self, "Expected exactly one alignments file.")

                root, ext = os.path.splitext(os.path.basename(input_paths[0]))
                # Complain if necessary
                if not ext =='.bam':
                    raise StepError(self, "The file %s does not appear to be any "
                                 "of bam.gz, bam.gzip, or bam"
                                 % input_paths[0]
                    )

                bedgraph_file = run.add_output_file(
                    'bedgraph', '%s.bg' % run_id, input_paths)
                bigwig_file = run.add_output_file(
                    'bigwig', '%s.bw' %run_id, input_paths)
                # Start creation of BedGraph files
                with run.new_exec_group() as bedgraph_group:
                    with bedgraph_group.add_pipeline() as pipe:
                        # BAM -> BedGraph
                        # (necessary for bedGraph, bigWig)
                        genomecov = [
                            self.get_tool('bedtools'), 'genomecov']
                        genomecov.append('-bg')
                        genomecov.append('-ibam')
                        genomecov.extend(input_paths )

                        pipe.add_command(genomecov)

                        sort = [ self.get_tool('sort'), '-k1,1', '-k2,2n']
                        pipe.add_command(
                            sort, stdout_path = bedgraph_file)
                with run.new_exec_group() as bigwig_group:
                    bedgraph_to_bigwig = [
                        self.get_tool('bedGraphToBigWig'),
                        bedgraph_file,
                        self.get_option('chromosome-sizes'),
                        bigwig_file ]
                    bigwig_group.add_command(bedgraph_to_bigwig)
