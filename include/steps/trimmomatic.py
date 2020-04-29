import os
from abstract_step import AbstractStep


class Trimmomatic(AbstractStep):
    def __init__(self, pipeline):
        super(Trimmomatic, self).__init__(pipeline)

        # Input
        # Single-end
        self.add_input_connection('first_read')

        # Paired-end
        self.add_input_connection('second_read')

        # Output
        # Single-end
        self.add_output_connection('log')
        self.add_output_connection('log_stderr')
        self.add_output_connection('forward')

        # Paired-end
        self.add_output_connection('forward.unpaired')
        self.add_output_connection('reverse')
        self.add_output_connection('reverse.unpaired')

        # Tools
        self.require_tool('java')

        # Options
        # Mandatory
        self.add_option('steps', list, optional=False,
                        description="List defining the analysis, in terms of "
                        "trimmomatic-commands")
        # Optional
        self.add_option('threads', int, optional=True, default=1,
                        description="Number of threads to use")
        self.add_option('phred-base', int, optional=True, choices=[33, 64],
                        description="PHRED-base")
        self.add_option('base_file_name', int, optional=True,
                        description="The prefix commoon to all output files, "
                        "replacing the run-ID")
        self.add_option('jar_path', str, optional=True, default="trimmomatic",
                        description="Path to trimmomatic.jar")

    def runs(self, run_ids_connections_files):
        for run_id in run_ids_connections_files.keys():
            trimmomatic_base = [self.get_tool('java')]
            trimmomatic_base.extend(('-jar', self.get_option('jar_path')))

            if run_ids_connections_files[run_id]['in/second_read'][0]:
                trimmomatic_base.append('PE')
                paired_end = True
            else:
                trimmomatic_base.append('SE')
                paired_end = False

            if self.is_option_set_in_config('threads'):
                trimmomatic_base.extend(('-threads',
                                         str(self.get_option('threads'))))

            if self.is_option_set_in_config('phred-base'):
                if self.get_option('phred-base') == 64:
                    trimmomatic_base.append('-phred64')
            else:
                trimmomatic_base.append('-phred33')

            step_list = list()
            for option in self.get_option('steps'):
                option[0] = option[0].upper()
                if option[0] == 'ILLUMINACLIP':
                    assert(len(option) in range(4, 7))
                    if option[1][-6:] != '.fasta':
                        raise Exception("Adapters should be provided as "
                                        "FASTA-file.")
                    if int(option[2]) not in list(range(17)):
                        raise Exception("The alignment-seed can contain up "
                                        "to 17 mismatches.")
                elif option[0] == 'SLIDINGWINDOW':
                    assert(len(option[1:]) == 2)
                    if int(option[2]) not in list(range(33, 158)):
                        raise Exception("SLIDINGWINDOW: Specify a valid "
                                        "PHRED-score.")
                elif option[0] == 'MAXINFO':
                    assert(len(option[1:]) == 2)
                    if not 0 <= float(option[2]) <= 1:
                        raise Exception("Strictness should be set between "
                                        " 0.0 and 1.0")
                elif option[0] == 'LEADING':
                    assert(len(option[1:]) == 1)
                    if int(option[1]) not in list(range(33, 158)):
                        raise Exception("LEADING: Specify a valid "
                                        "PHRED-score.")
                elif option[0] == 'TRAILING':
                    assert(len(option[1:]) == 1)
                    if int(option[1]) not in list(range(33, 158)):
                        raise Exception("TRAILING: Specify a valid "
                                        "PHRED-score.")
                elif option[0] == 'CROP':
                    assert(len(option[1:]) == 1)
                    assert(str(option[1]).isdigit())
                elif option[0] == 'HEADCROP':
                    assert(len(option[1:]) == 1)
                    assert(str(option[1]).isdigit())
                elif option[0] == 'MINLEN':
                    assert(len(option[1:]) == 1)
                    assert(str(option[1]).isdigit())
                elif option[0] == 'AVGQUAL':
                    assert(len(option[1:]) == 1)
                    if int(option[1]) not in list(range(33, 158)):
                        raise Exception("Specify a valid PHRED-score.")
                elif option[0] == 'TOPHRED33':
                    assert(len(option) == 1)
                elif option[0] == 'TOPHRED64':
                    assert(len(option) == 1)
                else:
                    raise Exception("TRIMMOMATIC: Invalid step declared")
                for i in range(len(option)):
                    option[i] = str(option[i])
                step_list.append(':'.join(option))

            with self.declare_run(run_id) as run:
                if paired_end:
                    if len(run_ids_connections_files[run_id]
                           ['in/first_read']) !=\
                        len(run_ids_connections_files[run_id]
                            ['in/second_read']):
                        raise Exception("Incorrect pairing of paired-end-"
                                        "files in run %s" % run_id)
                for file_no in range(len(run_ids_connections_files[run_id]
                                         ['in/first_read'])):

                    if self.is_option_set_in_config('base_file_name'):
                        file_id = self.get_option('base_file_name')
                    else:
                        file_id = run_id
                    if len(run_ids_connections_files[run_id]
                           ['in/first_read']) > 1:
                        file_id += str(file_no + 1)

                    with run.new_exec_group() as trimmomatic_exec_group:
                        option_list = list()
                        option_list.extend(['-trimlog', '%s.log' % file_id])

                        option_list.append(
                            run_ids_connections_files[run_id]
                            ['in/first_read'][file_no])
                        if paired_end:
                            option_list.append(
                                run_ids_connections_files[run_id]
                                ['in/second_read'][file_no])

                        option_list.append('%s.forward.fastq' % file_id)

                        if paired_end:
                            pe_out_files = ('forward.unpaired',
                                            'reverse',
                                            'reverse.unpaired')
                            for output_file in pe_out_files:
                                option_list.append('%s.%s.fastq' %
                                                   (file_id, output_file))

                        trimmomatic = trimmomatic_base[:]
                        trimmomatic.extend(option_list)
                        trimmomatic.extend(step_list)
                        trimmomatic_exec_group.add_command(
                            trimmomatic,
                            stderr_path=run.add_output_file(
                                'log_stderr',
                                '%s.log_stderr.txt' % file_id,
                                [run_ids_connections_files[run_id]
                                 ['in/first_read'][file_no]]))

                        run.add_output_file(
                            'log',
                            '%s.log' % file_id,
                            [run_ids_connections_files[run_id]
                             ['in/first_read'][file_no]])

                        run.add_output_file(
                            'forward',
                            '%s.forward.fastq' % file_id,
                            [run_ids_connections_files[run_id]
                             ['in/first_read'][file_no]])

                        if paired_end:
                            for output_file in pe_out_files:
                                run.add_output_file(
                                    output_file,
                                    '%s.%s.fastq' % (file_id, output_file),
                                    [run_ids_connections_files[run_id]
                                     ['in/first_read'][file_no]] +
                                    [run_ids_connections_files[run_id]
                                     ['in/second_read'][file_no]])
