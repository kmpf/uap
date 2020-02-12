import sys
from abstract_step import *
import glob
import misc
import process_pool
import yaml
import os

from logging import getLogger

logger=getLogger('uap_logger')

class gffCompare(AbstractStep):
    
    '''
    gffcompare [-r <reference_mrna.gtf> [-R]] [-G] [-T] [-V] [-s <seq_path>]
        [-o <outprefix>] [-p <cprefix>]
        {-i <input_gtf_list> | <input1.gtf> [<input2.gtf> .. <inputN.gtf>]}

     GffCompare provides classification and reference annotation mapping and
     matching statistics for RNA-Seq assemblies (transfrags) or other generic
     GFF/GTF files.
     GffCompare also clusters and tracks transcripts across multiple GFF/GTF
     files (samples), writing matching transcripts (identical intron chains) into
     <outprefix>.tracking, and a GTF file <outprefix>.combined.gtf which
     contains a nonredundant set of transcripts across all input files (with
     a single representative transfrag chosen for each clique of matching transfrags
     across samples).
    '''


    def __init__(self, pipeline):
        super(gffCompare, self).__init__(pipeline)

        self.set_cores(2)

        self.add_connection('in/assembling')
        self.add_connection('out/combined')
        self.add_connection('out/loci')  # *.loci
        self.add_connection('out/stats')      # *.stats
        self.add_connection('out/tracking')   # *.tracking        
        self.add_connection('out/log_stderr')

        self.require_tool('mkdir')
        self.require_tool('mv')
        self.require_tool('gffcompare')

        self.add_option('i', str, optional=True,
                        description='rrovide a text file with a list of (query) GTF files to process instead of expecting them as command line arguments useful when a large number of GTF files should be processed)')
        self.add_option('r', str, optional=True,
                        description='reference annotation file (GTF/GFF)')
        self.add_option('R', bool, optional=True,
                        description='for -r option, consider only the reference transcripts that overlap any of the input transfrags (Sn correction)')
        self.add_option('Q', bool, optional=True,
                        description='for -r option, consider only the input transcripts that overlap any of the reference transcripts (Precision correction); (Warning: this will discard all "novel" loci!)')
        self.add_option('M', bool, optional=True,
                        description='discard (ignore) single-exon transfrags and reference transcripts')
        self.add_option('N', bool, optional=True,
                        description='discard (ignore) single-exon reference transcripts')
        self.add_option('s', str, optional=True,
                        description='path to genome sequences (optional); this can be either a multi-FASTA file or a directory containing single-fasta files (one for each contig); repeats must be soft-masked (lower case) in order to be able to classify transfrags as repeats')
        self.add_option('e', int, optional=True,
                        description='max. distance (range) allowed from free ends of terminal exons of reference transcripts when assessing exon accuracy (100)')
        self.add_option('d', int, optional=True,
                        description='max. distance (range) for grouping transcript start sites (100)')
        self.add_option('p', str, optional=True,
                        description='rthe name prefix to use for consensus transcripts in the <outprefix>.combined.gtf file (default: \'TCONS\')')
        self.add_option('C', bool, optional=True,
                        description='include the "contained" transcripts in the .combined.gtf file')
        self.add_option('F', bool, optional=True,
                        description='do not discard intron-redundant transfrags if they share the 5\' end (if they differ only at the 3\' end)')
        self.add_option('G', bool, optional=True,
                        description='generic GFF input file(s): do not assume Cufflinks/Stringtie GTF input, (do not discard intron-redundant transfrags)')        
        self.add_option('T', bool, optional=True,
                        description='do not generate .tmap and .refmap files for each input file')  


    def runs(self, run_ids_connections_files):

        options=['i', 'r','R','Q','M','N','s','e','d','p','C','F','G','T']
        file_options=['i', 'r', 's']

        set_options = [option for option in options if \
                       self.is_option_set_in_config(option)]

        option_list = list()
        for option in set_options:
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    option_list.append('-%s' % option)
            else:
                value = str(self.get_option(option))
                if option in file_options:
                    value = os.path.abspath(value)
                option_list.append('-%s' % option)
                option_list.append(value)


        for run_id in run_ids_connections_files.keys():

             with self.declare_run(run_id) as run:

                base = run.get_output_directory_du_jour_placeholder() + '/' + run_id

                assembling  = run_ids_connections_files[run_id]['in/assembling'][0]
                input_paths = assembling

                combined = run.add_output_file(
                    'combined',
                    '%s.annotated.gtf' % run_id,
                    [input_paths])

                log_stderr = run.add_output_file(
                    'log_stderr',
                    '%s-log_stderr.txt' % run_id,
                    [input_paths])

                loci_file = run.add_output_file(
                    'loci',
                    '%s.loci' % run_id,
                    [input_paths])

                stats_file = run.add_output_file(
                    'stats',
                    '%s.stats' % run_id,
                    [input_paths])

                tracking_file = run.add_output_file(
                    'tracking',
                    '%s.tracking' % run_id,
                    [input_paths])

                # check, if only a single input file is provided
                len_input = run_ids_connections_files[run_id]['in/assembling']
                if len(len_input) != 1:
                    raise StandardError("Expected exactly one assembling file., but got this %s" % input_paths)

                with run.new_exec_group() as exec_group:

                   # exec_group.add_command(loci, stdout_path = loci_file)
                    with exec_group.add_pipeline() as pipe:
                        gffCompare = [self.get_tool('gffcompare'), '-o' , base]

                        gffCompare.extend(option_list)
                        gffCompare.append(assembling)

                        pipe.add_command(gffCompare, stdout_path=combined, stderr_path=log_stderr)

