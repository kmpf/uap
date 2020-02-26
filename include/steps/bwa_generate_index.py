from uaperrors import StepError
import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class BwaGenerateIndex(AbstractStep):
    '''
    This step generates the index database from sequences in the FASTA format.

    Typical command line::

        bwa index -p <index-basename> <seqeunce.fasta>
    '''

    def __init__(self, pipeline):
        super(BwaGenerateIndex, self).__init__(pipeline)
        self.set_cores(6)

        self.add_connection('in/reference_sequence')
        self.add_connection('out/bwa_index')

        self.require_tool('bwa')

        self.add_option('index-basename', str, optional = False,
                        description="Prefix of the created index database")

    def runs(self, run_ids_connections_files):

        for run_id in run_ids_connections_files.keys():
            # Get the basename
            index_basename =  "%s-%s" % (
                self.get_option('index-basename'), run_id)

            with self.declare_run(index_basename) as run:
                with run.new_exec_group() as exec_group:
                    refseq = run_ids_connections_files[run_id]\
                             ['in/reference_sequence']

                    if refseq == [None]:
                        raise StepError(self, "No reference sequence received.")
                    if len(refseq) != 1:
                        raise StepError(self,
                            "Reference sequence is not a single file.")
                    bwa_index = [self.get_tool('bwa'), 'index']
                    # Add index_basename
                    bwa_index.extend(
                        ['-p', os.path.join(
                            run.get_output_directory_du_jour_placeholder(),
                            index_basename)]
                    )
                    # Add reference sequence (a single file)
                    bwa_index.append(refseq[0])
                    exec_group.add_command(bwa_index)

                    run.add_output_file(
                        'bwa_index',
                        '%s.amb' % index_basename,
                        refseq)
                    run.add_output_file(
                        'bwa_index',
                        '%s.ann' % index_basename,
                        refseq)
                    run.add_output_file(
                        'bwa_index',
                        '%s.bwt' % index_basename,
                        refseq)
                    run.add_output_file(
                        'bwa_index',
                        '%s.pac' % index_basename,
                        refseq)
                    run.add_output_file(
                        'bwa_index',
                        '%s.sa' % index_basename,
                        refseq)
