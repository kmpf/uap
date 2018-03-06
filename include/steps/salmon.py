import os
from logging import getLogger
from abstract_step import AbstractStep

logger = getLogger('uap_logger')


class Salmon(AbstractStep):
    '''
    Salmon


    '''

    def __init__(self, pipeline):
        super(Salmon, self).__init__(pipeline)

        # input connections
        self.add_connection('in/first_read')
        self.add_connection('in/second_read')

        # output connections
        self.add_connection('out/cmd_info.json')
        self.add_connection('out/lib_format_counts.json')
        self.add_connection('out/quant.sf')
        self.add_connection('out/quant.genes.sf')
        self.add_connection('out/log_stderr')
        self.add_connection('out/log_stdout')

        self.dir_files = {
            'flenDist.txt': 'libParams',
            'salmon_quant.log': 'logs',
            'ambig_info.tsv': 'aux_info',
            'expected_bias.gz': 'aux_info',
            'fld.gz': 'aux_info',
            'meta_info.json': 'aux_info',
            'observed_bias.gz': 'aux_info',
            'observed_bias_3p.gz': 'aux_info'
        }

        for dir_file in self.dir_files:
            self.add_connection('out/' + dir_file)

        # required tools
        self.require_tool('salmon')
        self.require_tool('mkdir')
        self.require_tool('mv')
        self.require_tool('rm')

        # options
        self.add_option('cores', int, optional=True, default=1,
                        description="workaround to specify cores for grid \
                        engine and threads ie")

        self.add_option('i', str, optional=False, default=None,
                        description="Salmon index")

        self.add_option('g', str, optional=True, default=None,
                        description="File containing a mapping of transcripts " \
                                    "to genes.  If this file is provided " \
                                    "Salmon will output both quant.sf and " \
                                    "quant.genes.sf files, where the latter " \
                                    "contains aggregated gene-level abundance " \
                                    "estimates. The transcript to gene mapping " \
                                    "should be provided as either a GTF file, " \
                                    "or a in a simple tab-delimited format " \
                                    "where each line contains the name of a " \
                                    "transcript and the gene to which it " \
                                    "belongs separated by a tab. The " \
                                    "extension of the file is used to " \
                                    "determine how the file should be parsed. " \
                                    "Files ending in '.gtf', '.gff' or '.gff3'" \
                                    "are assumed to be in GTF format; files " \
                                    "with any other extension are assumed to " \
                                    "be in the simple format. In GTF / GFF " \
                                    "format, the 'transcript_id' is assumed " \
                                    "to contain the transcript identifier and " \
                                    "the 'gene_id' is assumed to contain the " \
                                    "corresponding gene identifier.")

    def runs(self, run_ids_connections_files):
        self.set_cores(self.get_option('cores'))

        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                input_fileset = []
                r1 = run_ids_connections_files[run_id]['in/first_read'][0]
                input_fileset.append(r1)

                r2 = None
                if 'in/second_read' in run_ids_connections_files[run_id]:
                    r2 = run_ids_connections_files[run_id]['in/second_read'][0]
                    input_fileset.append(r2)

                # create tmp directory for output
                with run.new_exec_group() as salmon_eg:
                    temp_dir = run.add_temporary_directory('salmon-out')
                    mkdir = [self.get_tool('mkdir'), temp_dir]
                    salmon_eg.add_command(mkdir)

                    salmon = [self.get_tool('salmon'), 'quant']
                    salmon.extend(['-i', str(self.get_option('i'))])
                    salmon.extend(['-l', 'ISF'])

                    if self.is_option_set_in_config('g'):
                        salmon.extend(['-g', self.get_option('g')])

                    salmon.extend(['-o', temp_dir])

                    # TODO: implement param --auxDir

                    salmon.extend(['-1', r1])
                    (r2 is not None) and (salmon.extend(['-2', r2]))

                    stderr_file = "%s-salmon-log_stderr.txt" % (run_id)
                    log_stderr = run.add_output_file("log_stderr",
                                                     stderr_file, input_fileset)
                    stdout_file = "%s-salmon-log_stdout.txt" % (run_id)
                    log_stdout = run.add_output_file("log_stdout",
                                                     stdout_file, input_fileset)

                salmon_eg.add_command(salmon, stdout_path=log_stdout,
                                            stderr_path=log_stderr)

                result_files = dict()
                result_files["cmd_info.json"] = run.add_output_file(
                    "cmd_info.json", "cmd_info.json", input_fileset)
                result_files["lib_format_counts.json"] = run.add_output_file(
                    "lib_format_counts.json", "lib_format_counts.json", input_fileset)
                result_files["quant.sf"] = run.add_output_file(
                    "quant.sf", "%s_quant.sf" % run_id, input_fileset)

                if self.is_option_set_in_config('g'):
                    result_files["quant.genes.sf"] = run.add_output_file(
                        "quant.genes.sf", "%s_quant.genes.sf" % run_id, input_fileset)
                else:
                    run.add_empty_output_connection("quant.genes.sf")

                result_files["flenDist.txt"] = run.add_output_file(
                    "flenDist.txt", "flenDist.txt", input_fileset)
                result_files["salmon_quant.log"] = run.add_output_file(
                    "salmon_quant.log", "salmon_quant.log", input_fileset)
                result_files["ambig_info.tsv"] = run.add_output_file(
                    "ambig_info.tsv", "ambig_info.tsv", input_fileset)
                result_files["expected_bias.gz"] = run.add_output_file(
                    "expected_bias.gz", "expected_bias.gz", input_fileset)
                result_files["fld.gz"] = run.add_output_file(
                    "fld.gz", "afld.gz", input_fileset)
                result_files["meta_info.json"] = run.add_output_file(
                    "meta_info.json", "meta_info.json", input_fileset)
                result_files["observed_bias.gz"] = run.add_output_file(
                    "observed_bias.gz", "observed_bias.gz", input_fileset)
                result_files["observed_bias_3p.gz"] = run.add_output_file(
                    "observed_bias_3p.gz", "observed_bias_3p.gz", input_fileset)

                # move file from temp directory to expected position
                with run.new_exec_group() as mv_exec_group:
                    for orig, dest_path in result_files.iteritems():
                        orig_path = os.path.join(temp_dir, orig)
                        if orig in self.dir_files:
                            orig_path = os.path.join(temp_dir, self.dir_files[orig], orig)
                            # TODO: how to modify dest_path?
                        mv = [self.get_tool('mv'), orig_path, dest_path]
                        mv_exec_group.add_command(mv)

                # remove directories in temp
                with run.new_exec_group() as rm_exec_group:
                    dirs = [temp_dir + '/aux_info', temp_dir + '/libParams', temp_dir + '/logs']
                    rm = [self.get_tool('rm'), '-r']
                    rm.extend(dirs)
                    rm_exec_group.add_command(rm)