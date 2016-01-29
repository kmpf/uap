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
        self.add_output_connection('forward')
        self.add_output_connection('log')
        self.add_output_connection('log_stdout')

        # Paired-end
        self.add_output_connection('reverse')
        self.add_output_connection('forward.unpaired')
        self.add_output_connection('reverse.unpaired')

        # Tools
        self.require_tool('trimmomatic')

        # Options
        # Mandatory
        self.add_option('paired-end', bool, optional=False, default=True,
                        description="Single-end- or paired-end-mode")
        self.add_option('steps', list, optional=False, default={},
                        description="List defining the analysis")
        # Optional
        self.add_option('threads', int, optional=True, default=1,
                        description="Number of threads to use")

        self.add_option('phred-base', int, optional=True, choices=[33, 64],
                        description="PHRED-base")

    def runs(self, run_ids_connections_files):
        trimmomatic = [self.get_tool('trimmomatic')]

        if self.get_option('paired_end'):
            trimmomatic.append('PE')
        else:
            trimmomatic.append('SE')

        if self.get_option('phred-base') == 64:
            trimmomatic.append('phred64')
        else:
            trimmomatic.append('phred33')

        for run_id in run_ids_connections_files.keys():
            trimmomatic = [self.get_tool('trimmomatic')]
            step_list = list()
            for option in self.get_option('steps'):
                option = option.upper()
                if option[0] is 'ILLUMINACLIP':
                    assert(len(option[1:]) == 4)
                    assert(option[2].isdigit())
                    assert(option[3].isdigit())
                    assert(option[4].isdigit())
                elif option[0] is 'SLIDINGWINDOW':
                    assert(len(option[1:]) == 2)
                    assert(option[1].isdigit())
                    assert(option[2].isdigit())
                elif option[0] is 'MAXINFO':
                    assert(len(option) == 1)
                elif option[0] is 'LEADING':
                    assert(len(option[1:]) == 1)
                    assert(option[1].isdigit())
                elif option[0] is 'TRAILING':
                    assert(len(option[1:]) == 1)
                    assert(option[1].isdigit())
                elif option[0] is 'CROP':
                    assert(len(option[1:]) == 1)
                    assert(option[1].isdigit())
                elif option[0] is 'HEADCROP':
                    assert(len(option[1:]) == 1)
                    assert(option[1].isdigit())
                elif option[0] is 'MINLEN':
                    assert(len(option[1:]) == 1)
                    assert(option[1].isdigit())
                elif option[0] is 'AVGQUAL':
                    assert(len(option[1:]) == 1)
                    assert(option[1].isdigit())
                elif option[0] is 'TOPHRED33':
                    assert(len(option) == 1)
                elif option[0] is 'TOPHRED64':
                    assert(len(option) == 1)
                else:
                    raise StandardError("TRIMMOMATIC: Invalid step declared")
                step_list.append(':'.join(option))

            with self.declare_run(run_id) as run:
                if self.get_option('paired-end'):
                    if len(run_ids_connections_files[run_id]
                           ['in/first_read']) !=\
                        len(run_ids_connections_files[run_id]
                            ['in/second_read']):
                        raise StandardError("Incorrect pairing of paired-end-"
                                            "files in run %s" % run_id)
                for file_no in range(len(run_ids_connections_files[run_id]
                                         ['in/first_read'])):

                    if len(run_ids_connections_files[run_id]
                           ['in/first_read']) > 1:
                        file_id = '%s_%d' % (run_id, file_no + 1)
                    else:
                        file_id = run_id

                    with run.new_exec_group() as trimmomatic_exec_group:
                        option_list = list()
                        option_list.append('-trimlog')
                        option_list.append(
                            os.path.join(
                                run.get_output_directory_du_jour_placeholder(),
                                '%s.log' % file_id))

                        option_list.append(
                            run_ids_connections_files[run_id]
                            ['in/first_read'][file_no])
                        if self.get_option('paired_end'):
                            option_list.append(
                                run_ids_connections_files[run_id]
                                ['in/second_read'][file_no])

                        option_list.append(
                            os.path.join(
                                run.get_output_directory_du_jour_placeholder(),
                                '%s.forward.fastq' % file_id))
                        if self.get_option('paired_end'):
                            pe_out_files = ('forward.unpaired',
                                            'reverse',
                                            'reverse.unpaired')

                            for output_file in pe_out_files:
                                option_list.append(
                                    os.path.join(
                                     run.
                                     get_output_directory_du_jour_placeholder(),
                                     '%s.%s.fastq' % (file_id, output_file)))

                        option_list.extend(step_list)
                        trimmomatic.extend(option_list)
                        trimmomatic_exec_group.add_command(
                            trimmomatic,
                            stdout_path=run.add_output_file(
                                'stdout',
                                '%s.log_stdout.txt' % file_id,
                                run_ids_connections_files[run_id]
                                ['ini/first_read']))

                        output_files = dict()

                        output_files['log'] =\
                            run.add_output_file(
                                output_file,
                                '%s.log' % file_id,
                                run_ids_connections_files[run_id]
                                    ['in/first_read'][file_no])

                        output_files['forward'] =\
                            run.add_output_file(
                                output_file,
                                '%s.forward.fastq' % file_id,
                                run_ids_connections_files[run_id]
                                    ['in/first_read'][file_no])

                        if self.get_option('paired_end'):
                            pe_out_files = ('forward.unpaired',
                                            'reverse',
                                            'reverse.unpaired')
                            for output_file in pe_out_files:
                                output_files[output_file] =\
                                    run.add_output_file(
                                        output_file,
                                        '%s.%s.fastq' % (file_id, output_file),
                                        run_ids_connections_files[run_id]
                                        ['in/first_read'][file_no] +
                                        run_ids_connections_files[run_id]
                                        ['in/second_read'][file_no])
