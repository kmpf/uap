from uaperrors import UAPError
import sys
import os
from logging import getLogger
from abstract_step import AbstractStep

logger=getLogger('uap_logger')

class BwaMem(AbstractStep):
    '''
    Align 70bp-1Mbp query sequences with the BWA-MEM algorithm. Briefly, the
    algorithm works by seeding alignments with maximal exact matches (MEMs) and
    then extending seeds with the affine-gap Smith-Waterman algorithm (SW).

    http://bio-bwa.sourceforge.net/bwa.shtml

    Typical command line::

        bwa mem [options] <bwa-index> <first-read.fastq> [<second-read.fastq>] \
        > <sam-output>
    '''

    def __init__(self, pipeline):
        super(BwaMem, self).__init__(pipeline)
        self.set_cores(6)

        self.add_connection('in/first_read')
        self.add_connection('in/second_read')
        self.add_connection('out/alignments')


        self.require_tool('dd')
        self.require_tool('mkfifo')
        self.require_tool('pigz')
        self.require_tool('bwa')

        # Options to set bwa mem flags
        self.add_option('index', str, optional=False,
                        description="Path to BWA index")
        ## [Algorithm options:]
        self.add_option('t', int, optional = True, default = 6,
                        description = "number of threads [6]")
        self.add_option('k', int, optional = True,
                        description = "minimum seed length [19]")
        self.add_option('w', int, optional = True,
                        description = "band width for banded alignment [100]")
        self.add_option('d', int, optional = True,
                        description = "off-diagonal X-dropoff [100]")
        self.add_option('r', float, optional = True,
                        description = "look for internal seeds inside a seed "
                        "longer than {-k} * FLOAT [1.5]")
        self.add_option('y', int, optional = True,
                        description = "seed occurrence for the 3rd round "
                        "seeding [20]")
        self.add_option('c', int, optional = True,
                        description = "skip seeds with more than INT "
                        "occurrences [500]")
        self.add_option('D', float, optional = True,
                        description = "drop chains shorter than FLOAT fraction "
                        "of the longest overlapping chain [0.50]")
        self.add_option('W', int, optional = True,
                        description = "discard a chain if seeded bases shorter "
                        "than INT [0]")
        self.add_option('m', int, optional = True,
                        description = "perform at most INT rounds of mate "
                        "rescues for each read [50]")
        self.add_option('S', bool, optional = True,
                        description = "skip mate rescue")
        self.add_option('P', bool, optional = True,
                        description = "skip pairing; mate rescue performed "
                        "unless -S also in use")
        self.add_option('e', bool, optional = True,
                        description = "discard full-length exact matches")

        ## [Scoring options:]
        self.add_option('A', int, optional = True,
                        description = "score for a sequence match, which "
                        "scales options -TdBOELU unless overridden [1]")
        self.add_option('B', int, optional = True,
                        description = "penalty for a mismatch [4]")
        self.add_option('O', str, optional = True,
                        description = "gap open penalties for deletions and "
                        "insertions [6,6]")
        self.add_option('E', str, optional = True,
                        description = "gap extension penalty; a gap of size k "
                        "cost '{-O} + {-E}*k' [1,1]")
        self.add_option('L', str, optional = True,
                        description = "penalty for 5'- and 3'-end clipping "
                        "[5,5]")
        self.add_option('U', int, optional = True,
                        description = "penalty for an unpaired read pair [17]")
        self.add_option('x', str, optional = True,
                        description = "read type. Setting -x changes multiple "
                        "parameters unless overriden [null]::\n"
                        "       pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio "
                        "reads to ref)\n"
                        "       ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford "
                        "Nanopore 2D-reads to ref)\n"
                        "       intractg: -B9 -O16 -L5  (intra-species contigs to ref)"
        )

        ## [Input/output options:]
        self.add_option('p', bool, optional = True,
                        description = "smart pairing (ignoring in2.fq)")
        self.add_option('R', str, optional = True,
                        description = "read group header line such as "
                        "'@RG\tID:foo\tSM:bar' [null]")
        self.add_option('H', str, optional = True,
                        description = "insert STR to header if it starts with "
                        "@; or insert lines in FILE [null]")
        self.add_option('j', bool, optional = True,
                        description = "treat ALT contigs as part of the "
                        "primary assembly (i.e. ignore <idxbase>.alt file)")
        self.add_option('v', int, optional = True,
                        description = "verbose level: 1=error, 2=warning, "
                        "3=message, 4+=debugging [3]")
        self.add_option('T', int, optional = True,
                        description = "minimum score to output [30]")
        self.add_option('h', str, optional = True,
                        description = "if there are <INT hits with score >80% "
                        "of the max score, output all in XA [5,200]")
        self.add_option('a', bool, optional = True,
                        description = "output all alignments for SE or "
                        "unpaired PE")
        self.add_option('C', bool, optional = True,
                        description = "append FASTA/FASTQ comment to SAM output")
        self.add_option('V', bool, optional = True,
                        description = "output the reference FASTA header in "
                        "the XR tag")
        self.add_option('Y', str, optional = True,
                        description = "use soft clipping for supplementary "
                        "alignments")
        self.add_option('M', str, optional = True,
                        description = "mark shorter split hits as secondary")

        # Options for dd
        self.add_option('dd-blocksize', str, optional = True, default = "256k")

    def runs(self, run_ids_connections_files):

        # Check if index is valid
        if not os.path.exists(self.get_option('index') + '.bwt'):
            raise UAPError("Could not find index: %s.*" %
                         self.get_option('index') )

        # Compile the list of options
        options = [
            # [Algorithm options:]
            't', 'k', 'w', 'd', 'r', 'y', 'c', 'D', 'W', 'm', 'S', 'P', 'e',
            # [Scoring options:]
            'A', 'B', 'O', 'E', 'L', 'U', 'x',
            # [Input/output options:]
            'p', 'R', 'H', 'j', 'v', 'T', 'h', 'a', 'C', 'V', 'Y', 'M'
        ]

        set_options = [option for option in options if \
                       self.is_option_set_in_config(option)]

        option_list = list()
        for option in set_options:
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    option_list.append('-%s' % option)
            else:
                option_list.append('-%s' % option)
                option_list.append(str(self.get_option(option)))

        for run_id in run_ids_connections_files.keys():
            with self.declare_run(run_id) as run:
                # Get list of files for first/second read
                fr_input = run_ids_connections_files[run_id]['in/first_read']
                sr_input = run_ids_connections_files[run_id]['in/second_read']

                input_paths = [ y for x in [fr_input, sr_input] \
                               for y in x if y !=None ]

                # Do we have paired end data and is it exactly one ?
                is_paired_end = False if sr_input == [None] else True

                # Fail if we have don't have exactly one file or
                # an empty connection
                if len(fr_input) != 1 or fr_input == [None]:
                    raise UAPError("Expected single input file for first read.")
                # Fail if we don't have exactly one file
                if is_paired_end and len(sr_input) != 1:
                    raise UAPError("Expected single input file for second read.")
                input_paths = fr_input # single element list
                if is_paired_end:
                    input_paths.extend(sr_input)

                # Check file endings for proper type
                for input_path in input_paths:
                    if len([_ for _ in ['fastq', 'fq', 'fq.gz', 'fastq.gz']\
                               if input_path.endswith(_)]) != 1:
                        raise UAPError("%s possess unknown suffix. "
                                     "(None of: fastq, fq, fq.gz, fastq.gz)")
                # BWA can handle only single files for first and second read
                # IMPORTANT: BWA handles gzipped as well as not gzipped files

                with run.new_exec_group() as exec_group:

                    def prepare_input(input_path, exec_group):
                        # Create temporary fifo
                        temp_fifo = run.add_temporary_file(
                            'in-fifo-%s' %
                            os.path.basename(input_path) )
                        mkfifo = [self.get_tool('mkfifo'), temp_fifo]
                        exec_group.add_command(mkfifo)
                        dd = [self.get_tool('dd'),
                              'bs=%s' % self.get_option('dd-blocksize'),
                              'if=%s' % input_path,
                              'of=%s' % temp_fifo]
                        exec_group.add_command(dd)

                        return (exec_group, temp_fifo)

                    # Temporary fifos
                    temp_fr_fifo, temp_sr_fifo = (str, str)

                    exec_group, temp_fr_fifo = prepare_input(
                            fr_input[0], exec_group)
                    # And if we handle paired end data
                    if is_paired_end:
                        exec_group, temp_sr_fifo = prepare_input(
                            sr_input[0], exec_group)

                    # 3. Map reads using bwa mem
                    with exec_group.add_pipeline() as bwa_mem_pipe:
                        # Assemble bwa mem command
                        bwa_mem = [
                            self.get_tool('bwa'),
                            'mem',
                        ]
                        bwa_mem.extend(option_list)
                        bwa_mem.append(self.get_option('index'))
                        bwa_mem.append(temp_fr_fifo)

                        if is_paired_end:
                            bwa_mem.append(temp_sr_fifo)

                        bwa_mem_pipe.add_command(bwa_mem)
                        # Compress bwa mem output
                        pigz = [self.get_tool('pigz'),
                                '--stdout']
                        bwa_mem_pipe.add_command(pigz)
                        # Write bowtie2 output to file
                        dd = [
                            self.get_tool('dd'),
                            'obs=%s' % self.get_option('dd-blocksize'),
                            'of=%s' %
                            run.add_output_file(
                                'alignments',
                                '%s-bwa-mem.sam.gz' % run_id,
                                input_paths
                            )
                        ]
                        bwa_mem_pipe.add_command(dd)
