import sys
import os
import re
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class Cutadapt(AbstractStep):
    '''
    Cutadapt finds and removes adapter sequences, primers, poly-A tails and
    other types of unwanted sequence from your high-throughput sequencing reads.

    https://cutadapt.readthedocs.org/en/stable/


    '''
    
    def __init__(self, pipeline):
        super(Cutadapt, self).__init__(pipeline)
        
        self.set_cores(4)
        
        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/first_read')
        self.add_connection('out/second_read')
        self.add_connection('out/log_first_read')
        self.add_connection('out/log_second_read')
        
        self.require_tool('cat')
        self.require_tool('cutadapt')
        self.require_tool('dd')
        self.require_tool('fix_qnames')
        self.require_tool('mkfifo')
        self.require_tool('pigz')

        # Options for cutadapt
        self.add_option('adapter-type', str, optional = True, default='-a',
                        choices=['-a', '-g', '-b', ''],
                        description="")
        self.add_option('adapter-R1', str, optional = True,
                        description="Adapter sequence to be clipped off of the"
                        "first read.")
        self.add_option('adapter-R2', str, optional = True,
                        description="Adapter sequence to be clipped off of the"
                        "second read")
        self.add_option('adapter-file', str, optional = True,
                        description="File containing adapter sequences to be "
                        "clipped off of the reads.")
        self.add_option('use_reverse_complement', bool, default = False,
                        description="The reverse complement of adapter "
                        "sequences 'adapter-R1' and 'adapter-R2' are used for "
                        "adapter clipping.")
        self.add_option('fix_qnames', bool, default = False,
                        description="If set to true, only the leftmost string "
                        "without spaces of the QNAME field of the FASTQ data is "
                        "kept. This might be necessary for downstream analysis.")
        self.add_option('m', int, default = 1)
        self.add_option('q', int, default = None, optional =True)
        self.add_option('u', int, default = None, optional =True)

    def runs(self, run_ids_connections_files):

        ## Make sure the adapter type is one of -a, -b or -g 
        if self.is_option_set_in_config('adapter-type'):
            if not self.get_option('adapter-type') in set(['-a','-b','-g', '']):
                logger.error("Option 'adapter-type' must be either '-a', "
                             "'-b', or '-g'!")
                sys.exit(1)

        read_types = {'first_read': 'R1', 'second_read': 'R2'}
        paired_end_info = dict()
        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                for read in read_types:                
                    connection = 'in/%s' % read
                    input_paths = run_ids_connections_files[run_id][connection]

                    if input_paths == [None]:
                        run.add_empty_output_connection("%s" % read)
                        run.add_empty_output_connection("log_%s" % read)
                    else:
                        paired_end_info[run_id] = self.find_upstream_info_for_input_paths(input_paths, 'paired_end')
                        # make sure that adapter-R1/adapter-R2 or adapter-file are
                        # correctly set 
                        # this kind of mutual exclusive option checking is a bit 
                        # tedious, so we do it here.
                        if paired_end_info[run_id]:
                            if ( not self.is_option_set_in_config('adapter-R2') and 
                                 not self.is_option_set_in_config('adapter-file') ):
                                logger.error(
                                    "Option 'adapter-R2' or 'adapter-file' "
                                    "required because sample %s is paired end!"
                                    % run_id)
                                sys.exit(1)
                        elif ( self.is_option_set_in_config('adapter-R2') and
                               not self.is_option_set_in_config('adapter-file') ):
                            logger.error(
                                "Option 'adapter-R2' not allowed because "
                                "sample %s is not paired end!" % run_id)
                            sys.exit(1)
                        if ( self.is_option_set_in_config('adapter-file') and
                             self.is_option_set_in_config('adapter-R1') ):
                            logger.error(
                                "Option 'adapter-R1' and 'adapter-file' "
                                "are both set but are mutually exclusive!")
                            sys.exit(1)
                        if ( not self.is_option_set_in_config('adapter-file') and
                             not self.is_option_set_in_config('adapter-R1') ):
                            logger.error(
                                "Option 'adapter-R1' or 'adapter-file' "
                                "required to call cutadapt for sample %s!"
                                % run_id)
                            sys.exit(1)
                        temp_fifos = list()
                        exec_group = run.new_exec_group()
                        for input_path in input_paths:
                            # 1. Create temporary fifo for every input file
                            temp_fifo = run.add_temporary_file(
                                "fifo-%s" % os.path.basename(input_path) )
                            temp_fifos.append(temp_fifo)
                            mkfifo = [self.get_tool('mkfifo'), temp_fifo]
                            exec_group.add_command(mkfifo)
                            # 2. Output files to fifo
                            if input_path.endswith('fastq.gz'):
                                with exec_group.add_pipeline() as pigz_pipe:
                                    # 2.1 command: Read file in 4MB chunks
                                    dd_in = [self.get_tool('dd'),
                                           'ibs=4M',
                                           'if=%s' % input_path]
                                    # 2.2 command: Uncompress file to fifo
                                    pigz = [self.get_tool('pigz'),
                                            '--decompress',
                                            '--stdout']
                                    # 2.3 command: Write file in 4MB chunks to 
                                    #              fifo
                                    dd_out = [self.get_tool('dd'),
                                              'obs=4M',
                                              'of=%s' % temp_fifo]

                                    pigz_pipe.add_command(dd_in)
                                    pigz_pipe.add_command(pigz)
                                    pigz_pipe.add_command(dd_out)

                            elif input_path.endswith('fastq'):
                                # 2.1 command: Read file in 4MB chunks and
                                #              write to fifo in 4MB chunks
                                dd_in = [self.get_tool('dd'),
                                         'bs=4M',
                                         'if=%s' % input_path,
                                         'of=%s' % temp_fifo]
                                exec_group.add_command(dd_in)
                            else:
                                logger.error("File %s does not end with any "
                                             "expected suffix (fastq.gz or "
                                             "fastq). Please fix that issue.")
                                sys.exit(1)
                        # 3. Read data from fifos
                        with exec_group.add_pipeline() as cutadapt_pipe:
                            # 3.1 command: Read from ALL fifos
                            cat = [self.get_tool('cat')]
                            cat.extend(temp_fifos)
                            cutadapt_pipe.add_command(cat)
                            # 3.2 command: Fix qnames if user wants us to
                            if self.get_option('fix_qnames') == True:
                                fix_qnames = [self.get_tool('fix_qnames')]
                                cutadapt_pipe.add_command(fix_qnames)

                            # Let's get the correct adapter sequences or
                            # adapter sequence fasta file 
                            adapter = None
                            # Do we have adapter sequences as input?
                            if self.is_option_set_in_config('adapter-%s' \
                                                            % read_types[read]):
                                # Get adapter sequence 
                                adapter = self.get_option(
                                    'adapter-%s' % read_types[read])
                            
                                # add index to adapter sequence if necessary
                                if '((INDEX))' in adapter:
                                    index = self.find_upstream_info_for_input_paths(
                                        input_paths,
                                        'index-%s' % read_types[read])
                                    adapter = adapter.replace('((INDEX))', index)
                            
                                # create reverse complement if necessary
                                if self.get_option('use_reverse_complement'):
                                    complements = string.maketrans('acgtACGT',
                                                                   'tgcaTGCA')
                                    adapter = adapter.translate(complements)[::-1]
                                
                                # make sure the adapter is looking good
#                                if re.search('^[ACGT]+$', adapter) == None:
#                                    logger.error("Unable to come up with a "
#                                                 "legit-looking adapter: %s"
#                                                 % adapter)
#                                    sys.exit(1)
                            # Or do we have a adapter sequence fasta file?
                            elif self.is_option_set_in_config('adapter-file'):
                                adapter = "file:" + self.get_option(
                                    'adapter-file')
                                if not os.path.exists(
                                        self.get_option('adapter-file')):
                                    logger.error(
                                        "File %s containing adapter sequences "
                                        "does not exist."
                                        % self.get_option('adapter-file'))
                                    sys.exit(1)


                            # 3.3 command: Clip adapters
                            cutadapt = [self.get_tool('cutadapt')]

                            if self.get_option('adapter-type') == '':
                                cutadapt.extend(['-'])
                            else:
                                cutadapt.extend([self.get_option('adapter-type'), 
                                                 adapter, '-'])


                            if self.is_option_set_in_config('m'):
                                cutadapt.extend(['-m',  str(self.get_option('m'))])
                            if self.is_option_set_in_config('q'):
                                cutadapt.extend(['-q',  str(self.get_option('q'))])
                            if self.is_option_set_in_config('u'):
                                cutadapt.extend(['-u',  str(self.get_option('u'))])


                            cutadapt_log_file = run.add_output_file(
                                    'log_%s' % read,
                                    '%s-cutadapt-%s-log.txt'
                                    % (run_id, read_types[read]),
                                    input_paths)
                            # 3.4 command: Compress output
                            pigz = [self.get_tool('pigz'),
                                     '--blocksize', '4096',
                                     '--processes', '1',
                                     '--stdout']
                            # 3.5 command: Write to output file in 4MB chunks
                            clipped_fastq_file = run.add_output_file(
                                "%s" % read,
                                "%s_%s.fastq.gz" %
                                (run_id, read_types[read]),
                                input_paths)

                            dd = [self.get_tool('dd'),
                                  'obs=4M',
                                  'of=%s' % clipped_fastq_file]


                            cutadapt_pipe.add_command(cutadapt,
                                                      stderr_path =\
                                                      cutadapt_log_file)
                            cutadapt_pipe.add_command(pigz)
                            cutadapt_pipe.add_command(dd)

    def reports(self, run_ids_connection_files):

        # Imports are done here to prevent them from causing errors while
        # analysis is done
        import matplotlib.pyplot as plt
        import numpy as np
        import csv

        table_header = ['length\tcount\texpect\tmax.err\terror counts']

        run_id_read_data = dict()

        for run_id in run_ids_connection_files:
            run_id_read_data[run_id] = dict()
            logs = ['log_first_read', 'log_second_read']
            for log in [l for l in logs \
                        if run_ids_connection_files[run_id][l] != [] ]:

                run_id_read_data[run_id][log] = dict()
                for log_file in run_ids_connection_files[run_id][log]:
                    base, ext = os.path.splitext( os.path.basename(log_file) )
                    tables = dict()
                    found_table = 0
                    record = False
                    with open(log_file, 'r') as f:
                        for line in f:
                            line = line.rstrip()
                            if line in table_header:
                                found_table += 1
                                tables[found_table] = list()
                                record = True
                                if record and line == '':
                                    record = False
                                    if record:
                                        tables[found_table].append(line)

                    # Plot a graph for every table found
                    data = dict()
                    for table_nr in tables.keys():
                        data[table_nr] = {
                            'length': list(),
                            'count': list(),
                            'expect': list(),
                            'max_error': list(),
                            'error_counts': list()
                        }
                        reader = csv.DictReader(tables[table_nr], delimiter='\t')
                        for row in reader:
                            data[table_nr]['length'].append( row['length'] )
                            data[table_nr]['count'].append( row['count'] )
                            data[table_nr]['expect'].append( row['expect'] )
                            data[table_nr]['max_error'].append( row['max.err'] )
                            data[table_nr]['error_counts'].append(
                                row['error counts'] )

                    run_id_read_data[run_id][log][log_file] = data

        # Plot the data in one immense figure
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for run_id, log in run_id_read_data.items():
            for log_file, data in run_id_read_data[run_id][log].items():
                for nr in data.keys():
                    ax.plot
                    print(data)
                    plt.plot(
                        data['length'],
                        data['count'],
                        'ro',
                        data['length'],
                        data['expect'],
                        'b--')

                    plt.xlim(0, float(data['length'][-1]) + 1)
                    plt.ylim(0.01, float(max(data['expect'][0], data['count'])) + 100 )
                    plt.yscale('log')
                    print('Trying to save report to %s.png' % base)
                    plt.savefig('%s-%s.png' % (base, table_nr))
                    plt.close()
