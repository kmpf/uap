from uaperrors import UAPError
import sys
from logging import getLogger
import os

from abstract_step import AbstractStep

logger = getLogger('uap_logger')


class MergeFastxFiles(AbstractStep):
    '''
    This step merges all .fastq/a(.gz) files belonging to a certain sample.
    First and second read files are merged separately. The output files are
    gzipped.
    '''

    def __init__(self, pipeline):
        super(MergeFastxFiles, self).__init__(pipeline)

        self.set_cores(4)  # muss auch in den Decorator

        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/first_read')
        self.add_connection('out/second_read')
        self.add_connection('out/report')

        self.require_tool('pigz')
        self.require_tool('echo')

    def _getFastFormat(self, fast_file, is_gzipped):

        required_file_extensions = [
            '.fastq', '.fq', 'fnq', '.fasta', '.fa', '.fna']

        example_file = os.path.basename(fast_file)
        format_index = -2 if is_gzipped else -1
        fast_format = '.' + example_file.split('.')[format_index]

        if fast_format not in required_file_extensions:
            raise UAPError("File %s does not end with any "
                         "expected suffix (%s). Please fix that issue." %
                         (fast_file, ' | '.join(required_file_extensions)))

        fast_char = fast_format[-1]
        return fast_char

    def _prepareFileList(self, run_id, file_list):

        head_line = ['run_id', 'first_read', 'second_read']
        list_data = ["\t".join(head_line) + "\n"]
        has_second_read = True if file_list['second_read'] != [None] else False

        for i, input_file in enumerate(file_list['first_read']):
            line = [run_id, os.path.basename(file_list['first_read'][i])]

            if has_second_read:
                line.append(os.path.basename(file_list['second_read'][i]))

            list_data.append("\t".join(line) + "\n")

        list_data[-1] = list_data[-1].rstrip("\n")
        return list_data

    def runs(self, run_ids_connections_files):

        read_types = {'first_read': '_R1', 'second_read': '_R2'}
        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                file_list = {'first_read': [], 'second_read': []}

                for read in read_types:
                    con = 'in/%s' % read
                    file_list[read] = run_ids_connections_files[run_id][con]

                echo_input = list(file_list['first_read'])
                if file_list['second_read'] != [None]:
                    echo_input.extend(file_list['second_read'])

                echo_out = run.add_output_file(
                    'report',
                    '%s.txt' % (run_id),
                    echo_input)
                prepared_list = self._prepareFileList(run_id, file_list)
                input_data = ''.join(prepared_list)
                echo_command = [self.get_tool('echo'), input_data]
                echo_exec_group = run.new_exec_group()
                echo_exec_group.add_command(echo_command, stdout_path=echo_out)

                for read in read_types:
                    connection = 'in/%s' % read
                    input_paths = run_ids_connections_files[run_id][connection]

                    if input_paths == [None]:
                        run.add_empty_output_connection("%s" % read)
                    else:
                        exec_group = run.new_exec_group()

                        # todo: outsource this as function
                        is_gzipped = False
                        for input_path in input_paths:
                            file_extension = os.path.splitext(input_path)[1]
                            is_gzipped = True if file_extension\
                                in ['.gz', '.gzip'] else False

                        # get fast-Format from inputfile for outputfile
                        fast_format = self._getFastFormat(input_paths[0],
                                                          is_gzipped)
                        p_out = run.add_output_file(
                            '%s' % read,
                            "%s%s.fast%s.gz" %
                            (run_id, read_types[read], fast_format),
                            input_paths)

                        if is_gzipped:
                            with exec_group.add_pipeline() as unzip_pipe:
                                pigz_input = [self.get_tool('pigz'),
                                              '--decompress', '--stdout']
                                pigz_input.extend(input_paths)

                                pigz_output = [self.get_tool('pigz'),
                                               '--stdout', '/dev/stdin']

                                unzip_pipe.add_command(pigz_input)
                                unzip_pipe.add_command(pigz_output,
                                                       stdout_path=p_out)
                        else:

                                pigz_output = [self.get_tool('pigz'),
                                               '--stdout']
                                pigz_output.extend(input_paths)
                                exec_group.add_command(pigz_output,
                                                       stdout_path=p_out)
