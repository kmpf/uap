import sys
from abstract_step import *
import glob
import misc
import process_pool
import yaml
import os

from logging import getLogger

logger=getLogger('uap_logger')

class CuffLinks(AbstractStep):

    '''
    CuffLinks is part of the 'Cufflinks suite of tools' for
    differential expr. analysis of RNA-Seq data and their
    visualisation. This step applies the cufflinks tool which
    assembles transcriptomes from RNA-Seq data and quantifies their
    expression and produces .gtf files with these annotations.
    For details on cufflinks we refer to the author's webpage:

    http://cole-trapnell-lab.github.io/cufflinks/

    '''
    def __init__(self, pipeline):
        super(CuffLinks, self).__init__(pipeline)

        self.set_cores(6)
        
        self.add_connection('in/alignments')
        self.add_connection('out/features')
        self.add_connection('out/skipped')
        self.add_connection('out/genes-fpkm')
        self.add_connection('out/isoforms_fpkm')
        self.add_connection('out/log_stderr')
        
        self.require_tool('mkdir')
        self.require_tool('mv')
        self.require_tool('cufflinks')

        ## options for cufflinks program
        # [General options:]
        # --out-put-dir
        self.add_option('num-threads', int, optional=True,
                        description='number of threads used during analysis')
        self.add_option('seed', int, optional=True, 
                        description='value of random number generator seed')
        self.add_option('GTF', bool, optional=True, 
                        description='quantitate against reference transcript annotations')
        self.add_option('GTF-guide', bool, optional=True, 
                        description='use reference transcript annotation to guide assembly')
        self.add_option('mask-file', str, optional=True, 
                        description='ignore all alignment within transcripts in this file')
        self.add_option('frag-bias-correct', str, optional=True, 
                        description='use bias correction - reference fasta required')
        self.add_option('multi-read-correct', bool, optional=True, 
                        description='use \'rescue method\' for multi-reads (more accurate)')
        self.add_option('library-type', str, optional=False,
                        choices=['ff-firststrand','ff-secondstrand','ff-unstranded','fr-firststrand',
                                 'fr-secondstrand','fr-unstranded','transfrags'], 
                        description='library prep used for input reads')
        self.add_option('library-norm-method', str, choices=['classic-fpkm'], optional=True,
                        description='Method used to normalize library sizes')

        # [Advanced Abundance Estimation Options:]
        self.add_option('frag-len-mean', int, optional=True,
                        description='average fragment length (unpaired reads only)')
        self.add_option('frag-len-std-dev', int, optional=True, 
                        description='fragment length std deviation (unpaired reads only)')
        self.add_option('max-mle-iterations', int, optional=True, 
                        description='maximum iterations allowed for MLE calculation')
        self.add_option('compatible-hits-norm', bool, optional=True, 
                        description='count hits compatible with reference RNAs only')
        self.add_option('total-hits-norm', bool, optional=True, 
                        description='count all hits for normalization')
        self.add_option('num-frag-count-draws', int, optional=True, 
                        description='Number of fragment generation samples')
        self.add_option('num-frag-assign-draws', int, optional=True, 
                        description='Number of fragment assignment samples per generation')
        self.add_option('max-frag-multihits', str, optional=True, 
                        description='Maximum number of alignments allowed per fragment')
        self.add_option('no-effective-length-correction', bool, optional=True, 
                        description='No effective length correction')
        self.add_option('no-length-correction', bool, optional=True, 
                        description='No length correction')

        # [Advanced Assembly Options:]
        self.add_option('label', str, optional=True, 
                        description='assembled transcripts have this ID prefix')
        self.add_option('min-isoform-fraction', float, optional=True, 
                        description='suppress transcripts below this abundance level')
        self.add_option('pre-mrna-fraction', float, optional=True, 
                        description='suppress intra-intronic transcripts below this level')
        self.add_option('max-intron-length', int, optional=True,
                        description='ignore alignments with gaps longer than this')
        self.add_option('junc-alpha', float, optional=True, 
                        description='alpha for junction binomial test filter')
        self.add_option('small-anchor-fraction', float, optional=True, 
                        description='percent read overhang taken as \'suspiciously small\'')
        self.add_option('min-frags-per-transfrag', int, optional=True, 
                        description='minimum number of fragments needed for new transfrags')
        self.add_option('overhang-tolerance', int, optional=True, 
                        description='number of terminal exon bp to tolerate in introns')
        self.add_option('max-bundle-length', int, optional=True, 
                        description='maximum genomic length allowed for a given bundle')
        self.add_option('max-bundle-frags', int, optional=True, 
                        description='maximum fragments allowed in a bundle before skipping')
        self.add_option('min-intron-length', int, optional=True, 
                        description='minimum intron size allowed in genome')
        self.add_option('trim-3-avgcov-thresh', int, optional=True, 
                        description='minimum avg coverage required to attempt 3\' trimming')
        self.add_option('trim-3-dropoff-frac', float, optional=True, 
                        description='fraction of avg coverage below which to trim 3\' end')
        self.add_option('max-multiread-fraction', float, optional=True, 
                        description='maximum fraction of allowed multireads per transcript')
        self.add_option('overlap-radius', int, optional=True, 
                        description='maximum gap size to fill between transfrags (in bp)')

        # [Advanced Reference Annotation Guided Assembly Options:]

        self.add_option('no-faux-reads', bool, optional=True, 
                        description='disable tiling by faux reads')
        self.add_option('3-overhang-tolerance', int, optional=True, 
                        description='overhang allowed on 3\' end when merging with reference')
        self.add_option('intron-overhang-tolerance', int, optional=True, 
                        description='overhang allowed inside reference intron when merging')

        # [Advanced Program Behavior Options:]
        self.add_option('verbose', bool, optional=True, 
                        description='log-friendly verbose processing (no progress bar)')
        # In uap we set this option to TRUE per default, since we only need this info within a log file
        # it suppresses the progress
        #self.add_option('quiet', bool, optional=True, 
        #                description='log-friendly quiet processing (no progress bar)')
        self.add_option('no-update-check', bool, optional=True, 
                        description='do not contact server to check for update availability')

    def runs(self, run_ids_connections_files):
        
        # Compile the list of options
        options=['num-threads', 'seed', 'GTF', 'GTF-guide', 'mask-file', 'frag-bias-correct', 'multi-read-correct', 'library-type', 
                 'library-norm-method', 'frag-len-mean', 'frag-len-std-dev', 'max-mle-iterations', 'compatible-hits-norm', 
                 'total-hits-norm', 'num-frag-count-draws', 'num-frag-assign-draws', 'max-frag-multihits', 'no-effective-length-correction', 
                 'no-length-correction', 'label', 'min-isoform-fraction', 'pre-mrna-fraction', 'max-intron-length', 'junc-alpha', 
                 'small-anchor-fraction', 'min-frags-per-transfrag', 'overhang-tolerance', 'max-bundle-length', 'max-bundle-frags', 
                 'min-intron-length', 'trim-3-avgcov-thresh', 'trim-3-dropoff-frac', 'max-multiread-fraction', 'overlap-radius', 
                 'no-faux-reads', '3-overhang-tolerance', 'intron-overhang-tolerance', 'verbose', 'no-update-check']
#                 'no-faux-reads', '3-overhang-tolerance', 'intron-overhang-tolerance', 'verbose', 'quiet', 'no-update-check']

        set_options = [option for option in options if \
                       self.is_option_set_in_config(option)]

        option_list = list()
        for option in set_options:
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    option_list.append('--%s' % option)
            else:
                option_list.append( '--%s' % option )
                option_list.append( str(self.get_option(option)) )


        for run_id in run_ids_connections_files.keys():

             with self.declare_run(run_id) as run:
                input_paths = run_ids_connections_files[run_id]['in/alignments']
                temp_dir = run.add_temporary_directory('cufflinks-out')
                
                # check, if only a single input file is provided
                if len(input_paths) != 1:
                    raise StandardError("Expected exactly one alignments file., but got this %s" % input_paths)

#                # check if temporary directory is there 
#                if not os.path.isdir(tmp_dir):
#                    #dir not present
#                    logger.error("Directory %s not found" % tmp_dir)
#                    sys.exit(1)
#
#                # .. and if its accessible
#                if not os.access(tmp_dir, os.W_OK):
#                    #not accessible
#                    logger.error("Directory %s not accessible." % tmp_dir)
#                    sys.exit(1)

                cufflinks = [self.get_tool('cufflinks'),'-o', temp_dir, '-q']
                cufflinks.extend(option_list)
                cufflinks.append(input_paths[0])

                # 3. Announce output files
                #    !!!Keep in mind that run.add_output_file() returns the
                #    expected path of the output file
                # Das hier sind die Output files die du UAP bekannt geben willst
                # habsch geklaut aus altem Step ;)

                result_files = {
                    'transcripts.gtf': run.add_output_file(
                        'features',
                        '%s-transcripts.gtf' % run_id,
                        input_paths
                        ),
                    'skipped.gtf': run.add_output_file(
                        'skipped',
                        '%s-skipped.gtf' % run_id,
                        input_paths
                        ),
                    'genes.fpkm_tracking': run.add_output_file(
                        'genes-fpkm',
                        '%s-genes.fpkm_tracking' % run_id,
                        input_paths
                        ),
                    'isoforms.fpkm_tracking': run.add_output_file(
                        'isoforms_fpkm',
                        '%s-isoforms.fpkm_tracking' % run_id,
                        input_paths
                        )
                    }
             
             with run.new_exec_group() as exec_group:
                 # 1. Create temporary directory for cufflinks output
                 mkdir = [self.get_tool('mkdir'), temp_dir]
                 exec_group.add_command(mkdir)

             # 2. Execute cufflinks
             # So wie ich die Sache sehe braucht man kein add_pipeline
             #with exec_group.add_pipeline() as cufflinks_pipe:
             #    cufflinks_pipe.add_command(cufflinks)
             # Du willst ja nur das ein Kommando ausfuehren
             with run.new_exec_group() as exec_group:
                 exec_group.add_command(cufflinks,
                                        stderr_path = run.add_output_file('log_stderr',
                                                                          '%s-cufflinks-log.txt' % run_id,
                                                                          input_paths
                                                                          )#,
                                        #hints = {'writes': result_files.values()}
                                        )

             # Files get moved to expected position after cufflinks finished
             with run.new_exec_group() as mv_exec_group:
                 for orig, dest_path in result_files.iteritems():
                     # 3. Rename files 
                     orig_path = os.path.join(temp_dir, orig)
                     mv = [self.get_tool('mv'), orig_path, dest_path]
                     mv_exec_group.add_command(mv)
