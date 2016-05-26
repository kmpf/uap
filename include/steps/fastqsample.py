import logging
from abstract_step import AbstractStep
import os
import sys

logger = logging.getLogger('uap_logger')


class FastqSample(AbstractStep):
    '''
    wrapper class for fastq-sample
    sample random reads from a fastq file
    http://homes.cs.washington.edu/~dcjones/fastq-tools/fastq-sample.html

    for a specific seed the subsampling process will ever produce the
    same order of positions so the connections between R1 and R2 remains
    (paired end)
    '''

    def __init__(self, pipeline):
        super(FastqSample, self).__init__(pipeline)

        self.set_cores(1)  # muss auch in den Decorator

        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/first_read')
        self.add_connection('out/second_read')

        self.require_tool('fastq-sample')
        self.require_tool('pigz')
        self.require_tool('mv')
        self.require_tool('rm')

        # Options for fastq-sample
        self.add_option('n', int, default=1000, optional=True,
                        description="The number of reads to sample and output")
        self.add_option('p', float, default=None, optional=True,
                        description="The number of reads to sample in terms "
                        "of the proportion of total reads. "
                        "If sampling with replacement, this number "
                        "may be greater than 1.0")
        self.add_option('o', str, default=None, optional=True,
                        description="The filename prefix to which output "
                        "should be written. If single-end data is being "
                        "sampled, the output file is [PREFIX].fastq, "
                        "and with paired-end, [PREFIX].1.fastq "
                        "and [PREFIX].2.fastq")
        self.add_option('r', str, default=None, optional=True,
                        description="Sample with replacement")
        self.add_option('c', str, default=None, optional=True,
                        description="Output reads not included in the random "
                        "sample to a file (or files) with the given prefix. "
                        "By default, these reads are not output.")
        self.add_option('s', str, default='1234', optional=True,
                        description="Seed the random number generator. "
                        "Using the same seed on the same data set will "
                        "produce the same random sample.")

        self.possible_options = ['n', 'p', 'o', 'r', 'c', 's']

    def runs(self, run_ids_connections_files):

        isset_n = self.is_option_set_in_config('n')
        isset_p = self.is_option_set_in_config('p')

        if isset_n and isset_p:
            logger.error("Option n AND p are set in config.yaml. "
                         "Only one is allowed.")
            sys.exit(1)

        config_options = self.get_options()

        read_types = {'first_read': '_R1', 'second_read': '_R2'}
        for run_id in run_ids_connections_files.keys():
            new_run_id = run_id
            # create new run id if option o isset
            if self.is_option_set_in_config('o'):
               new_run_id = config_options['o'] + '_' + run_id

            with self.declare_run(new_run_id) as run:

                for read in read_types:
                    connection = 'in/%s' % read
                    input_paths = run_ids_connections_files[run_id][connection]

                    if input_paths == [None]:
                        run.add_empty_output_connection("second_read")
                    else:
                        for input_path in input_paths:
                            # Get base name of input file
                            root, ext = os.path.splitext(
                                os.path.basename(input_path))

                            temp_file = input_path

                            is_gzipped = False
                            file_ext = os.path.splitext(input_path)[1]
                            is_gzipped = True if file_ext\
                                in ['.gz', '.gzip'] else False

                            if is_gzipped:
                                parts = os.path.basename(input_path).split('.')
                                root = '.'.join(parts[:-2])

                                # Unzip fastq
                                temp_file = run.add_temporary_file()
                                pigz_decompress_eg = run.new_exec_group()
                                pigz = [self.get_tool('pigz'),
                                        '--decompress', '--keep',
                                        '--stdout', input_path]

                                pigz_decompress_eg.add_command(
                                    pigz, stdout_path=temp_file)
                            # 1. Run fastqc for input file
                            fastqsample_eg = run.new_exec_group()

                            # @todo: its impossible to get a shorter line at
                            # this position for pep8-compatibility...
                            # maybe rename method?
                            outfile_path = run.get_output_directory_du_jour_placeholder()
                            outfile = outfile_path + "/sample"

                            fastqsample = [self.get_tool('fastq-sample')]

                            for option, value in config_options.iteritems():
                                if option in self.possible_options:
                                    if option == 'o':
                                        continue
                                    fastqsample.extend(['-%s' % (option),
                                                       str(value)])

                            fastqsample.extend(['-o', outfile])
                            fastqsample.append(temp_file)
                            fastqsample_eg.add_command(fastqsample)

                            #output compress subsample
                            filename_params = (new_run_id,
                                               read_types[read])

                            subsample_file = run.add_output_file(
                                "%s" % read,
                                "%s%s.fastq.gz" % filename_params,
                                [input_path])

                            pigz_compress_eg = run.new_exec_group()
                            pigz_compress = [self.get_tool('pigz'),
                                             '--best', '--stdout',
                                             outfile + '.fastq']
                            pigz_compress_eg.add_command(
                                pigz_compress, stdout_path=subsample_file)

                            # deletions
                            remove_eg = run.new_exec_group()
                            remove = [self.get_tool('rm'),
                                      outfile + '.fastq']
                            remove_eg.add_command(remove)
