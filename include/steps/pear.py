import os
from abstract_step import AbstractStep


class Pear(AbstractStep):
    def __init__(self, pipeline):
        super(Pear, self).__init__(pipeline)

        self.add_input_connection('first_read')
        self.add_input_connection('second_read')

        self.add_output_connection('assembled')
        self.add_output_connection('unassembled.forward')
        self.add_output_connection('unassembled.reverse')
        self.add_output_connection('discarded')
        self.add_output_connection('log')

        self.require_tool('pear')

        self.add_option('p-value', float, optional=True, default=0.01,
                        choices=[0.0001, 0.001, 0.01, 0.05, 1.0],
                        description='Specify a p-value for the statistical '
                        'test.')

        self.add_option('min-overlap', int, optional=True, default=10,
                        description='Specify the minimum overlap size.')

        self.add_option('m', int, optional=True,
                        default=300, description='Specify the maximum '
                        'possible length of the assembled sequences.')

        self.add_option('n', int, optional=True,
                        default=50, description='Specify the minimum '
                        'possible length of the assembled sequences.')

        self.add_option('min-trim-length', int, optional=True, default=1,
                        description='Specify the minimum length of reads '
                        'after trimming the low quality part.')

        self.add_option('q', int, optional=True, default=0,
                        choices=[0] + list(range(33, 158)),
                        description='Specify the quality score threshold '
                        'for trimming the low quality part of a read.')

        self.add_option('max-uncalled-base', float, optional=True,
                        choices=['[0.0..1.0]'], default=1.0,
                        description='Specify the maximal proportion of '
                        'uncalled bases in a read.')

        self.add_option('test-method', int, optional=True, default=1,
                        choices=[1, 2], description='Specify the type of '
                        'statistical tesit.')

        self.add_option('empirical-freqs', bool, optional=True,
                        default=False, description='Disable empirical base'
                        'frequencies.')

        self.add_option('score-method', int, optional=True, default=2,
                        choices=[1, 2, 3], description='Specify the scoring '
                        'method.')

        self.add_option('phred-base', int, optional=True, default=33,
                        choices=[33, 64], description='Base PHRED quality '
                        'score.')

        self.add_option('memory', str, optional=True,
                        description='Specify the amount of memory to be used.')

        self.add_option('cap', int, optional=True, default=40,
                        choices=[0] + list(range(33, 158)),
                        description='Specify the upper bound for the '
                        'resulting quality score.')

        self.add_option('threads', int, optional=True,
                        description='Number of threads to use.')

        self.add_option('nbase', bool, optional=True, default=False,
                        description='When resolving a mismatching base-pair '
                        'out of which non is degenerate, set the merged base '
                        'to N and use the highest quality score of the two '
                        'bases.')

    def runs(self, run_ids_connections_files):
        param_list = self.get_param_list()
        out_files = ('assembled', 'unassembled.forward', 'unassembled.reverse',
                     'discarded')
        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                if len(run_ids_connections_files[run_id]['in/first_read']) !=\
                   len(run_ids_connections_files[run_id]['in/second_read']):
                    raise Exception("Incorrect pairing of paired-end-files "
                                        "in run %s" % run_id)
                for file_no in range(len(run_ids_connections_files[run_id]
                                         ['in/first_read'])):
                    option_list = list()
                    if len(run_ids_connections_files[run_id]
                           ['in/first_read']) > 1:
                        file_id = '%s_%d' % (run_id, file_no + 1)
                    else:
                        file_id = run_id

                    option_list.append('--forward-fastq')
                    option_list.append(run_ids_connections_files[run_id]
                                       ['in/first_read'][file_no])
                    option_list.append('--reverse-fastq')
                    option_list.append(run_ids_connections_files[run_id]
                                       ['in/second_read'][file_no])

                    with run.new_exec_group() as pear_exec_group:
                        option_list.append('--output')
                        option_list.append(file_id)

                        pear = [self.get_tool('pear')]
                        pear.extend(option_list)
                        pear.extend(param_list)
                        pear_exec_group.add_command(
                            pear,
                            stdout_path=run.add_output_file(
                                'log',
                                '%s.log_stdout.txt' % file_id,
                                run_ids_connections_files[run_id]
                                ['in/first_read'] +
                                run_ids_connections_files[run_id]
                                ['in/second_read']))

                        for output_file in out_files:
                            run.add_output_file(
                                    output_file,
                                    '%s.%s.fastq' % (file_id, output_file),
                                    [run_ids_connections_files[run_id]
                                    ['in/first_read'][file_no],
                                    run_ids_connections_files[run_id]
                                    ['in/second_read'][file_no]])
