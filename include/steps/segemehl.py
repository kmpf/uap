import sys

from abstract_step import AbstractStep
import misc

class Segemehl(AbstractStep):
    '''
    segemehl is a software to map short sequencer reads to reference genomes. 
    Unlike other methods, segemehl is able to detect not only mismatches but 
    also insertions and deletions. Furthermore, segemehl is not limited to a 
    specific read length and is able to mapprimer- or polyadenylation 
    contaminated reads correctly.

    This step creates at first two FIFOs. The first through which segemehl gets 
    its genome data and the second to which it writes unmapped reads::

       mkfifo genome_fifo unmapped_fifo
       cat <genome-fasta> -o genome_fifo

    The executed segemehl command is this::

        segemehl -d genome_fifo -i <genome-index-file> -q <read1-fastq> 
                 [-p <read2-fastq>] -u unmapped_fifo -H 1 -t 11 -s -S -D 0
                 -o /dev/stdout |  pigz --blocksize 4096 --processes 2 -c

    The unmapped reads are saved via these commands::

        cat unmapped_fifo | pigz --blocksize 4096 --processes 2 -c > 
        <unmapped-fastq>

    '''

    def __init__(self, pipeline):
        super(Segemehl, self).__init__(pipeline)

        self.set_cores(12)
        
        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/alignments')
        self.add_connection('out/unmapped')
        self.add_connection('out/log')
        
        self.require_tool('dd')
        self.require_tool('pigz')
        self.require_tool('segemehl')

        # Options to set segemehl flags
        ## [INPUT]
        self.add_option('genome', str)
        self.add_option('index', str)
        self.add_option('bisulfite', int, choices = [1, 2], optional = True)
        ## [GENERAL]
        self.add_option('minsize', int, optional = True)
        self.add_option('silent', bool, default = True, optional = True)
        self.add_option('brief', bool, default = False, optional = True)
        ## [SEEDPARAMS]
        self.add_option('differences', int, default = 1, optional = True)
        self.add_option('jump', int, optional = True)
        self.add_option('evalue', float, optional = True)
        self.add_option('maxsplitevalue', float, optional = True)
        self.add_option('maxinterval', int, optional = True)
        self.add_option('splits', bool, default = True, optional = True)
        self.add_option('SEGEMEHL', bool, optional = True)
        self.add_option('MEOP', optional = True)
        self.add_option('nohead', optional = True)
        ## [SEEDEXTENSIONPARAMS]
        self.add_option('extensionscore', int, optional = True)
        self.add_option('extensionpenalty', int, optional = True)
        self.add_option('dropoff', int, optional = True)

        ## [ALIGNPARAMS]
        self.add_option('accuracy', int, optional = True)
        self.add_option('minsplicecover', int, optional = True)
        self.add_option('minfragscore', int, optional = True)
        self.add_option('minfraglen', int, optional = True)
        self.add_option('splicescorescale', float, optional = True)
        self.add_option('hitstrategy', int, choices = [0, 1], optional = True,
                        default = 1)
        self.add_option('showalign', bool, optional = True)
        self.add_option('prime5', str, optional = True)
        self.add_option('prime3', str, optional = True)
        self.add_option('clipacc', int, optional = True)
        self.add_option('polyA', bool, optional = True)
        self.add_option('autoclip', bool, optional = True)
        self.add_option('hardclip', bool, optional = True)
        self.add_option('order', bool, optional = True)
        self.add_option('maxinsertsize', int, optional = True)

    def runs(self, run_ids_connections_files):
        # Compile the list of options
        options = ['bisulfite', 'minsize', 'silent', 'brief', 'differences',
                   'jump', 'evalue', 'maxsplitevalue', 'maxinterval', 'splits',
                   'SEGEMEHL', 'MEOP', 'nohead', 'extensionscore',
                   'extensionpenalty', 'dropoff', 'accuracy', 'minsplicecover',
                   'minfragscore', 'minfraglen', 'splicescorescale',
                   'hitstrategy', 'showalign', 'prime5', 'prime3', 'clipacc',
                   'polyA', 'autoclip', 'hardclip', 'order', 'maxinsertsize']

        set_options = [option for option in options if \
                       self.is_option_set_in_config(option)]

        option_list = list()
        for option in set_options:
            if self.get_option(option):
                option_list.append('--%s' % option)
            if not isinstance(self.get_option(option), bool):
                option_list.append(str(self.get_option(option)))

        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                # Get list of files for first/second read
                fr_input = run_ids_connections_files[run_id]['in/first_read']
                sr_input = run_ids_connections_files[run_id]['in/second_read']

                input_paths = [ y for x in [fr_input, sr_input] \
                               for y in x if y != None ]

                # Do we have paired end data?
                is_paired_end = True if sr_input == [None] else False

                first_read_file = 
                if len(run_ids_connections_files[run_id]['in/first_read']) == 1:
                    first_read_file = run_ids_connections_files[run_id]\
                                      ['in/first_read']
                else:
                    raise StandardError("Expected one input file.")

                if run_ids_connections_files[run_id]['in/second_read'] == [None]:
                    second_read_files = None
                elif len(run_ids_connections_files[run_id]['in/second_read']) == 1:
                    second_read_files = run_ids_connections_files[run_id]\
                                        ['in/second_read']
                else:
                    raise StandardError("Expected one input file.")

                if not os.path.isfile(self.get_option('index')):
                    raise StandardError(
                        "The path %s provided to 'index' is not a file."
                        % self.get_option('index') )

                # SEGEMEHL can cope with gzipped files so we do not need to!!!
                #is_fr_gzipped = True if os.path.splitext(first_read_file[0])[1]\
                #                 in ['.gz', '.gzip'] else False

                #is_sr_gzipped = True if os.path.splitext(second_read_file[0])[1]\
                #                 in ['.gz', '.gzip'] else False

                # Segemehl is run in this exec group
                # Can segemehl handle multiple input files/fifos?
                with run.new_exec_group() as exec_group:
                    # 1. Create FIFO for genome
                    fifo_path_genome = run.add_temporary_file(
                        'segemehl-genome-fifo', designation = 'input')
                    mkfifo_genome = [self.get_tool('mkfifo'), fifo_path_genome]
                    exec_group.add_command(mkfifo_genome)
                    # 2. Create FIFO to write unmapped reads to
                    fifo_path_unmapped = run.add_temporary_file(
                        'segemehl-unmapped-fifo', designation = 'output')
                    mkfifo_unmapped = [self.get_tool('mkfifo'), fifo_path_unmapped]
                    exec_group.add_command(mkfifo_unmapped)
                    # 3. Read genome and output to FIFO
                    dd_genome = [self.get_tool('dd'),
                                 'bs=4M',
                                 'if=%s' % self.get_option('genome')
                                 'of=%s' % fifo_path_genome]
                    exec_group.add_command(dd_genome)
                
                    with exec_group.new_pipeline() as segemhel_pipe:
                        # 4. Start segemehl
                        segemehl = [
                            self.get_tool('segemehl'),
                            '--database', fifo_path_genome,
                            '--index', self.get_option('index'),
                            '--nomatchfilename', fifo_path_unmapped,
                            '--threads', '10'
                            '--query', first_read_path
                        ]
                        if is_paired_end:
                            segemehl.extend(['--mate', second_read_path])
                        segemehl.extend(option_list)
                        segemehl_pipe.add_command(
                            segemehl,
                            stderr_path = run.add_output_file(
                                'log',
                                '%s-segemehl-log.txt' % run_id,
                                input_paths)
                        )

                        # 5. Compress segemehl mapped reads
                        pigz_mapped_reads = [
                            self.get_tool('pigz'),
                            '--stdout',
                            '--blocksize', '4096', 
                            '--processes', '1'
                        ]

                        segemehl_pipe.add_command(
                            pigz_mapped_reads,
                            stdout_path = run.add_output_file(
                                'alignments',
                                '%s-segemehl-results.sam.gz' % run_id,
                                input_paths)
                        )


                    # 6. Compress unmapped reads
                    pigz_unmapped_reads = [
                        self.get_tool('pigz'),
                        '--stdout',
                        '--blocksize', '4096',
                        '--processes', '1',
                        fifo_path_unmapped
                    ]
                    exec_group.add_command(
                        pigz_unmapped_reads,
                        stdout_path = run.add_output_file(
                            'unmapped',
                            '%s-segemehl-unmapped.fastq.gz' % run_id,
                            input_paths)
                    )
