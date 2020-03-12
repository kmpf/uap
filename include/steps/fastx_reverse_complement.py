from uaperrors import StepError
import logging
from abstract_step import AbstractStep
import os

logger = logging.getLogger('uap_logger')


class FastxReverseComplement(AbstractStep):
    '''
    wrapper class for fastx_reverse_complement from fastx toolkit
    creates reverse complement of fasta and fastq files.
    http://hannonlab.cshl.edu/fastx_toolkit/
    '''

    def __init__(self, pipeline):
        super(FastxReverseComplement, self).__init__(pipeline)

        self.set_cores(1)  # muss auch in den Decorator

        self.add_connection('in/fastx')
        self.add_connection('out/fastx')

        self.require_tool('fastx_reverse_complement')
        self.require_tool('pigz')
        self.require_tool('cat')

        self.add_option('prefix', str, default=None, optional=True,
                        description="Add Prefix to sample name (deprecated).")

        self.add_option('name sheme', str, optional=True,
                        default='%s_revcom',
                        description=r"Naming sheme for the output files "
                        r"without '.fastq.gz' extension and where ``%s`` "
                        r"is replaced with the run id.")

    def runs(self, run_ids_connections_files):
        run_id_sheme = self.get_option('name sheme')
        prefix = self.get_option('prefix')
        if prefix:
            run_id_sheme = '%s_%%s_R1' % prefix
            logger.warning("[%s] The 'prefix' option is deprecaded in favor "
                           "of the 'name sheme' option. The set pefix '%s' is "
                           "converted to 'name sheme: %s'" %
                           (self, prefix, run_id_sheme))
        try:
            _ = run_id_sheme % ''
        except TypeError as e:
            raise StepError(self, 'Could not parse name sheme "%s": %s' %
                            (run_id_sheme, e))
        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                input_paths = run_ids_connections_files[run_id]['in/fastx']
                if input_paths == [None]:
                    run.add_empty_output_connection("alignments")
                elif len(input_paths) != 1:
                    raise StepError(
                        self, "Expected exactly one alignments file.")
                else:
                    is_gzipped = True if os.path.splitext(input_paths[0])[1]\
                        in ['.gz', '.gzip'] else False

                out = run.add_output_file(
                    "fastx",
                    run_id_sheme % run_id + '.fastq.gz',
                    input_paths)

                with run.new_exec_group() as exec_group:
                    with exec_group.add_pipeline() as pipe:
                        # 1.1 command: Uncompress file
                        if is_gzipped:
                            pigz = [self.get_tool('pigz'),
                                    '--decompress',
                                    '--processes', '1',
                                    '--stdout']
                            pigz.extend(input_paths)
                            pipe.add_command(pigz)
                        else:
                            cat = [self.get_tool('cat')]
                            cat.extend(input_paths)
                            pipe.add_command(cat)

                        # 1. Run  fastx  for input file
                        fastx_revcom = [
                            self.get_tool('fastx_reverse_complement')]
                        # gzip
                        fastx_revcom.extend(['-z'])
                        pipe.add_command(fastx_revcom, stdout_path=out)
