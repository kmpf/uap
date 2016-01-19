from abstract_step import AbstractStep
import os


class Pear(AbstractStep):
    def __init__(self, pipeline):
        super(Pear, self).__init__(pipeline)

        self.add_input_connection('first_read')
        self.add_input_connection('second_read')

        self.add_output_connection('read')
        self.add_output_connection('first_read')
        self.add_output_connection('second_read')
        self.add_output_connection('discarded_read')

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
                        choices=[0] + range(33, 158),
                        description='Specify the quality score threshold '
                        'for trimming the low quality part of a read.')

        self.add_option('max-uncalled-base', float, optional=True,
                        choices=['[0.0..1.0]'],
                        default=1.0, description='Specify the maximal '
                        'proportion of uncalled bases in a read.')

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
                        choices=[0] + range(33, 158),
                        description='Specify the upper bound for the resulting '
                        'quality score.')

        self.add_option('threads', int, optional=True,
                        description='Number of threads to use.')

        self.add_option('nbase', bool, optional=True, default=False,
                        description='When resolving a mismatching base-pair '
                        'out of which non is degenerate, set the merged base '
                        'to N and use the highest quality score of the two '
                        'bases.')

    def runs(self, run_ids_connections_files):
        param_list = list()
        for option, value in self.get_options('tool-specific').iteritems():
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    param_list.append('--' + str(option))
            else:
                if len(option) > 1:
                    param_list.append('--' + str(option))
                else:
                    param_list.append('-' + str(option))
                param_list.append(str(value))

        for run_id in run_ids_connections_files:
            with self.declare_run(run_id) as run:
                option_list = list()
                option_list.append('--forward-fastq')
                option_list.append(run_ids_connections_files[run_id]
                                   ['in/first_read'][0])
                option_list.append('--reverse-fastq')
                option_list.append(run_ids_connections_files[run_id]
                                   ['in/second_read'][0])

                with run.new_exec_group() as pear_exec_group:
                    option_list.append('--output')
                    option_list.append(
                        os.path.join(
                            run.get_output_directory_du_jour_placeholder(),
                            run_id))

                    pear = [self.get_tool('pear')]
                    option_list.extend(param_list)
                    pear.extend(option_list)
                    pear_exec_group.add_command(pear)

                    output_files = dict()
                    output_files['%s.assembled.fastq' % run_id] =\
                        run.add_output_file(
                            'read',
                            '%s.assembled.fastq' % run_id,
                            run_ids_connections_files[run_id]['in/first_read'] +
                            run_ids_connections_files[run_id]['in/second_read'])

                    output_files['%s.unassembled.forward.fastq' % run_id] =\
                        run.add_output_file(
                            'first_read',
                            '%s.unassembled.forward.fastq' % run_id,
                            run_ids_connections_files[run_id]['in/first_read'])

                    output_files['%s.unassembled.reverse.fastq' % run_id] =\
                        run.add_output_file(
                            'second_read',
                            '%s.unassembled.reverse.fastq' % run_id,
                            run_ids_connections_files[run_id]['in/second_read'])

                    output_files['%s.discarded.fastq' % run_id] =\
                        run.add_output_file(
                            'discarded_read',
                            '%s.discarded.fastq' % run_id,
                            run_ids_connections_files[run_id]['in/first_read'] +
                            run_ids_connections_files[run_id]['in/second_read'])
