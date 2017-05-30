import sys
from abstract_step import *
import glob
import misc
import process_pool
import yaml
import os

from logging import getLogger

logger=getLogger('uap_logger')

class CuffCompare(AbstractStep):

    '''
    CuffCompare is part of the 'Cufflinks suite of tools' for
    differential expr. analysis of RNA-Seq data and their
    visualisation. This step compares a cufflinks assembly to 
    known annotation. For details about cuffcompare we refer to
    the author's webpage:

    http://cole-trapnell-lab.github.io/cufflinks/cuffcompare/

    '''

    def __init__(self, pipeline):
        super(CuffCompare, self).__init__(pipeline)

        self.set_cores(1)

        self.add_connection('in/features')    # cuffmerge output
        self.add_connection('out/features')   # *.combined.gft
        self.add_connection('out/loci')       # *.loci
        self.add_connection('out/stats')      # *.stats
        self.add_connection('out/tracking')   # *.tracking
        self.add_connection('out/log_stderr')

        self.require_tool('cuffcompare')
#        self.require_tool('cd')

        # reference files
        self.add_option('r', str, optional=True,
                        description='An optional "reference" annotation GFF file '
                        'containing a set of known mRNAs to use as a reference '
                        'for assessing the accuracy of mRNAs or gene models '
                        'given in <input.gtf>') 
        self.add_option('s', str, optional=True,
                        description='Can be a multi-fasta file with all the genomic '
                        'sequences or a directory containing multiple single-fasta '
                        'files (one file per contig); lower case bases will be used '
                        'to classify input transcripts as repeats. NOTE that must '
                        'contain one fasta file per reference chromosome, and each '
                        'file must be named after the chromosome, and have a .fa or '
                        '.fasta extension.')

        # options
        self.add_option('R', bool, optional=True, # -R
                        description='For "-r" option, consider only '
                        'the reference transcripts that overlap any of the input '
                        'transfrags (Sn-correction)')
        self.add_option('Q', bool, optional=True, # -Q
                        description='For "-r" option, consider only '
                        'the input transcripts that overlap any of the reference '
                        'transcripts (Sp-correction)')
        self.add_option('M', bool, optional=True, # -M
                        description='Discard (ignore) single-exon transfrags and '
                        'reference transcripts')
        self.add_option('N', bool, optional=True, # -N
                        description='Discard (ignore) single-exon reference ')
        self.add_option('C', bool, optional=True, # -C
                        description='Enables the “contained” transcripts to be also '
                        'written in the .combined.gtffile, with the attribute '
                        '"contained_in" showing the first container transfrag found. '
                        'By default, without this option, cuffcompare does not write '
                        'in that file isoforms that were found to be fully '
                        'contained/covered (with the same compatible intron '
                        'structure) by other transfrags in the same locus.'
                        'transcripts')
        self.add_option('F', bool, optional=True, # -F
                        description='Do not discard intron-redundant transfrags if '
                        'they share the 5p end (if they differ only at the 3p end)')
        self.add_option('G', bool, optional=True,
                        description='generic GFF input file(s): do not assume '
                        'Cufflinks GTF, do not discard any intron-redundant '
                        'transfrags')
        self.add_option('V', bool, optional=True,
                        description='verbose processing mode (showing all GFF '
                        'parsing warnings)')
        
        # parameters
        self.add_option('e', int, optional=True, # -e
                        description='Max. distance (range) allowed from free ends '
                        'of terminal exons of reference transcripts when assessing '
                        'exon accuracy (Default: 100)')
        self.add_option('d', int, optional=True, # -d
                        description='Max. distance (range) for grouping transcript '
                        'start sites (Default: 100)')

        # skipped options are:
        # -T, do not generate .tmap and .refmap files for each input file
        #     --> uap expects these files and it is not allowed to suppress their
        #         generation

        # skipped paramters are:
        # -o <out-prefix> 
        # -c <consensus-prefix> 
        # --> both are set by uap and can not be defined by the usr!

    def runs(self, run_ids_connections_files):
        
        # Compile the list of options
        options = ['r', 's', 'R', 'Q', 'M', 'N', 'C', 'F', 'G', 'V', 'e', 'd']
            
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

                input_paths = run_ids_connections_files[run_id]['in/features']
                if not input_paths:
                    raise StandardError("No input files for run %s" % (run_id))
                    
                # check whether there's exactly one feature file
                if len(input_paths) != 1:
                    raise StandardError("Expected exactly one feature file.")

                in_file = input_paths[0]

                # the temporary output directory
                outdir = run.get_output_directory_du_jour_placeholder()

                # this is the prefix for the cufflinks cmd options:
                # -o and -c (out- and consensus prefix):
                prefixCC = '%s/%s_cuffcompare' % (outdir, run_id)

                # this is the prefix without directory for uap add output file
                prefix = '%s_cuffcompare' %  run_id
                            

                features_file = run.add_output_file('features',
                                                    '%s.combined.gtf' % prefix,
                                                    input_paths)
                loci_file     = run.add_output_file('loci',
                                                    '%s.loci' % prefix,
                                                    input_paths)
                stats_file    = run.add_output_file('stats',
                                                    '%s.stats' % prefix,
                                                    input_paths)
                tracking_file = run.add_output_file('tracking',
                                                    '%s.tracking' % prefix,
                                                    input_paths)
                log_err_file  = run.add_output_file('log_stderr',
                                                    '%s-cuffcompare-log_stderr.txt' % prefix,
                                                    input_paths)

                # create cuffcompare command
                # i) add fix options and parameters
                cuffcompare = [selfl.get_tool('cuffcompare'), 
                               '-o', prefixCC,
                               '-c', prefixCC,
                               '-T', # consider if we really want to suppress the map 
                                     # files or if we want to trace them
                               ]
                # ii) add user defined settings
                cuffcompare.extend(option_list)
                # iii) add input file
                cuffcompare.extend(in_file)

                with run.new_exec_group() as cc_exec_group:
                    cc_exec_group.add_command(cuffcompare)
