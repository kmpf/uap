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
        self.add_output_connection('log_stdout')
        self.add_output_connection('forward')

        # Paired-end
        self.add_output_connection('forward.unpaired')
        self.add_output_connection('reverse')
        self.add_output_connection('reverse.unpaired')

        # Tools
        self.require_tool('java')

        # Options
        # Mandatory
        self.add_option('jar_path', str, optional=False,
                        description="Path to trimmomatic.jar")
        self.add_option('paired-end', bool, optional=False,
                        description="Single-end- or paired-end-mode")
        self.add_option('steps', list, optional=False,
                        description="List defining the analysis, in terms of "
                        "trimmomatic-commands")
        # Optional
        self.add_option('threads', int, optional=True, default=1,
                        description="Number of threads to use")

        self.add_option('phred-base', int, optional=True, choices=[33, 64],
                        description="PHRED-base")

    def runs(self, run_ids_connections_files):
        trimmomatic_base = [self.get_tool('java')]
        trimmomatic_base.extend(('-jar', self.get_option('jar_path')))
        if self.get_option('paired-end'):
            trimmomatic_base.append('PE')
        else:
            trimmomatic_base.append('SE')

        if self.is_option_set_in_config('threads'):
            trimmomatic_base.extend(('-threads',
                                     str(self.get_option('threads'))))

        if self.is_option_set_in_config('phred-base'):
            if self.get_option('phred-base') == 64:
                trimmomatic_base.append('-phred64')
            else:
                trimmomatic_base.append('-phred33')

        for run_id in run_ids_connections_files.keys():
            step_list = list()
            for option in self.get_option('steps'):
                option[0] = option[0].upper()
                if option[0] == 'ILLUMINACLIP':
                    assert(len(option[1:]) == 4)
                    assert(str(option[2]).isdigit())
                    assert(str(option[3]).isdigit())
                    assert(str(option[4]).isdigit())
                elif option[0] == 'SLIDINGWINDOW':
                    assert(len(option[1:]) == 2)
                    assert(str(option[1]).isdigit())
                    assert(str(option[2]).isdigit())
                elif option[0] == 'MAXINFO':
                    assert(len(option) == 1)
                elif option[0] == 'LEADING':
                    assert(len(option[1:]) == 1)
                    assert(str(option[1]).isdigit())
                elif option[0] == 'TRAILING':
                    assert(len(option[1:]) == 1)
                    assert(str(option[1]).isdigit())
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
                    assert(str(option[1]).isdigit())
                elif option[0] == 'TOPHRED33':
                    assert(len(option) == 1)
                elif option[0] == 'TOPHRED64':
                    assert(len(option) == 1)
                else:
                    raise StandardError("TRIMMOMATIC: Invalid step declared")
                for i in range(len(option)):
                    option[i] = str(option[i])
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
                        file_id = '%s%d' % (run_id, file_no + 1)
                    else:
                        file_id = run_id

                    with run.new_exec_group() as trimmomatic_exec_group:
                        option_list = list()
                        option_list.extend(
                            ('-trimlog',
                             os.path.join(
                                 run.get_output_directory_du_jour_placeholder(),
                                 '%s.log' % file_id)))

                        option_list.append(
                            run_ids_connections_files[run_id]
                            ['in/first_read'][file_no])
                        if self.get_option('paired-end'):
                            option_list.append(
                                run_ids_connections_files[run_id]
                                ['in/second_read'][file_no])

                        option_list.append(
                            os.path.join(
                                run.
                                get_output_directory_du_jour_placeholder(),
                                '%s.forward.fastq' % file_id))

                        if self.get_option('paired-end'):
                            pe_out_files = ('forward.unpaired',
                                            'reverse',
                                            'reverse.unpaired')
                            for output_file in pe_out_files:
                                option_list.append(
                                    os.path.join(
                                     run.
                                     get_output_directory_du_jour_placeholder(),
                                     '%s.%s.fastq' % (file_id, output_file)))

                        trimmomatic = trimmomatic_base[:]
                        trimmomatic.extend(option_list)
                        trimmomatic.extend(step_list)
                        trimmomatic_exec_group.add_command(
                            trimmomatic,
                            stdout_path=run.add_output_file(
                                'log_stdout',
                                '%s.log_stdout.txt' % file_id,
                                [run_ids_connections_files[run_id]
                                 ['in/first_read'][file_no]]))

                        run.add_output_file(
                            output_file,
                            '%s.log' % file_id,
                            [run_ids_connections_files[run_id]
                             ['in/first_read'][file_no]])

                        run.add_output_file(
                            output_file,
                            '%s.forward.fastq' % file_id,
                            [run_ids_connections_files[run_id]
                             ['in/first_read'][file_no]])

                        if self.get_option('paired-end'):
                            for output_file in pe_out_files:
                                run.add_output_file(
                                    output_file,
                                    '%s.%s.fastq' % (file_id, output_file),
                                    [run_ids_connections_files[run_id]
                                     ['in/first_read'][file_no]] +
                                    [run_ids_connections_files[run_id]
                                     ['in/second_read'][file_no]])
