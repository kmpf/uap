import os
from logging import getLogger
from abstract_step import AbstractStep

logger = getLogger('uap_logger')

class ChromHmmLearnModel(AbstractStep):
    '''
    This command takes a directory with a set of binarized data files and learns
    a chromatin state model. Binarized data files have ‘_binary’ in the file
    name. The format for the binarized data files are that the first line
    contains the name of the cell separated by a tab with the name of the
    chromosome. The second line contains in tab delimited form the name of each
    mark. The remaining lines correspond to consecutive bins on the chromosome.
    The remaining lines in tab delimited form corresponding to each mark, with a
    ‘1’ for a present call or ‘0’ for an absent call and a ‘2’ if the data is
    considered missing at that interval for the mark.
    '''

    def __init__(self, pipeline):
        super(ChromHmmLearnModel, self).__init__(pipeline)

        self.set_cores(8)
        
        self.add_connection('in/cellmarkfiletable')
        self.add_connection('in/chromhmm_binarization')
        self.add_connection('out/')
        
        self.require_tool('ChromHMM')
        self.require_tool('echo')
        self.require_tool('ln')
        
        # ChromHMM LearnModel Required Parameters
        self.add_option('numstates', int, optional = False,
                        descritpion = "This parameter specifies the number of "
                        "states to include in the model.")
        self.add_option('assembly', str, optional = False,
                        description = "specifies the genome assembly. overlap "
                        "and neighborhood enrichments will be called with "
                        "default parameters using this genome assembly."
                        "Assembly names are e.g. hg18, hg19, GRCh38")
        # ChromHMM LearnModel Optional Parameters
        # [-b binsize]
        self.add_option('b', int, optional = True,
                        description = "The number of base pairs in a bin "
                        "determining the resolution of the model learning and "
                        "segmentation. By default this parameter value is set "
                        "to 200 base pairs.")
        # [-color r,g,b]
        self.add_option('color', str, optional = True,
                        description = 'This specifies the color of the heat '
                        'map. "r,g,b" are integer values between 0 and 255 '
                        'separated by commas. By default this parameter value '
                        'is 0,0,255 corresponding to blue.')
        # [-d convergedelta]
        self.add_option('d', float, optional = True,
                        description = "The threshold on the change on the "
                        "estimated log likelihood that if it falls below this "
                        "value, then parameter training will terminate. If "
                        "this value is less than 0 then it is not used as part "
                        "of the stopping criteria. The default value for this "
                        "parameter is 0.001.")
        # [-e loadsmoothemission]
        self.add_option('e', float, optional = True,
                        description = "This parameter is only applicable if "
                        "the load option is selected for the init parameter. "
                        "This parameter controls the smoothing away from 0 "
                        "when loading a model. The emission value used in the "
                        "model initialization is a weighted average of the "
                        "value in the file and a uniform probability over the "
                        "two possible emissions. The value in the file gets "
                        "weight (1-loadsmoothemission) while uniform gets "
                        "weight loadsmoothemission. The default value of this "
                        "parameter is 0.02.")
        # [-f inputfilelist]: Unnecessary due to the use of inputdir
        # [-h informationsmooth]
        self.add_option('h', float, optional = True,
                        description = "A smoothing constant away from 0 for "
                        "all parameters in the information based "
                        "initialization. This option is ignored if random or "
                        "load are selected for the initialization method. The "
                        "default value of this parameter is 0.02.")
        # [-holdcolumnorder]
        self.add_option('holdcolumnorder', bool, optional = True,
                        description = "Including this flag suppresses the "
                        "reordering of the mark columns in the emission "
                        "parameter table display.")
        # [-i outfileID]: If set is set programmatically
        # [-init information|random|load]
        self.add_option('init', str, optional = True,
                        choices = ["information","random","load"],
                        description = "This specifies the method for parameter "
                        "initialization method. 'information' is the default "
                        "method described in (Ernst and Kellis, Nature Methods "
                        "2012). 'random' – randomly initializes the parameters "
                        "from a uniform distribution. 'load' loads the "
                        "parameters specified in '-m modelinitialfile' and "
                        "smooths them based on the value of the "
                        "'loadsmoothemission' and 'loadsmoothtransition' "
                        "parameters. The default is information.")
        # [-l chromosomelengthfile] Should this be Mandatory???
        self.add_option('l', str, optional = True,
                        description = 'This file specifies the length of the '
                        'chromosomes. It is a two column tab delimited file '
                        'with the first column specifying the chromosome name '
                        'and the second column the length. If this file is '
                        'provided then no end coordinate will exceed what is '
                        'specified in this file. By default BinarizeBed '
                        'excludes the last partial bin along the chromosome, '
                        'but if that is included in the binarized data input '
                        'files then this file should be included to give a '
                        'valid end coordinate for the last interval.')
        # [-m modelinitialfile]
        self.add_option('m', str, optional = True,
                        description = "This specifies the model file "
                        "containing the initial parameters which can then be "
                        "used with the load option")
        # [-nobed]
        self.add_option('nobed', bool, optional = True,
                        description = "If this flag is present, then this "
                        "suppresses the printing of segmentation information "
                        "in the four column format. The default is to generate "
                        "a four column segmentation file")
        # [-nobrowser]
        self.add_option('nobrowser', bool, optional = True,
                        description = "If this flag is present, then browser "
                        "files are not printed. If -nobed is requested then "
                        "browserfile writing is also suppressed.")
        # [-noenrich]
        self.add_option('noenrich', bool, optional = True,
                        description = "If this flag is present, then "
                        "enrichment files are not printed. If -nobed is "
                        "requested then enrichment file writing is also "
                        "suppressed.")
        # [-p maxprocessors]
        self.add_option('p', int, optional = True,
                        description = "If this option is present then ChromHMM "
                        "will attempt to train a model using multiple "
                        "processors in parallel. ChromHMM will create multiple "
                        "threads up to the number of processors specified by "
                        "maxprocessors. If maxprocessors is set to 0, then the "
                        "number of threads will just be limited by the number "
                        "of processors available. Note that each additional "
                        "processor used will require additional memory. If "
                        "this option is specified ChromHMM will use a standard "
                        "Baum-Welch training algorithm opposed to the default "
                        "incremental expectation-maximization algorithm.")
#        # [-printposterior]
#        self.add_option('printposterior', bool, optional = True,
#                        description = "If this flag is present the posterior "
#                        "probabilities over state assignments are also printed "
#                        "in a file. These files end with ‘_posterior.txt’. One "
#                        "file is generated per cell type and chromosome. The "
#                        "first line of these files specify the chromosome and "
#                        "cell type, followed by a header line for each column, "
#                        "and then the posterior probabilities one per line. By "
#                        "default these files are not printed.")
#        # [-printstatebyline]
#        self.add_option('printstatebyline', bool, optional = True,
#                        description = "If this flag is present the state "
#                        "assignment are printed to a file one per line. These "
#                        "files end with ‘_maxstate.txt’. One file is generated "
#                        "per cell type and chromosome. The first line "
#                        "specifies the cell type and chromosome and the second "
#                        "line says MaxState and the state ordering methods. "
#                        "The remaining lines have the state assignments. By "
#                        "default these files are not printed.")

        # [-r maxiterations]
        self.add_option('r', int, optional = True,
                        description = "This option specifies the maximum "
                        "number of iterations over all the input data in the "
                        "training. By default this is set to 200.")
        # [-s seed]
        self.add_option('s', int, optional = True,
                        description = "This allows the specification of the "
                        "random seed. Randomization is used to determine the "
                        "visit order of chromosomes in the incremental "
                        "expectation-maximization algorithm used to train the "
                        "parameters and also used to generate the initial "
                        "values of the parameters if random is specified for "
                        "the init method.")
        # [-stateordering emission|transition]
        self.add_option('stateordering', str, optional = True,
                        choices = ["emission", "transition"],
                        description = "This determines whether the states are "
                        "ordered based on the emission or transition "
                        "parameters. See (Ernst and Kellis, Nature Methods) "
                        "for details. Default is 'emission'.")
        # [-t loadsmoothtransition]
        self.add_option('t', float, optional = True,
                        description = "This parameter is only applicable if "
                        "the load option is selected for the init parameter. "
                        "This parameter controls the smoothing away from 0 "
                        "when loading a model. The transition value used in "
                        "the model initialization is a weighted average of the "
                        "value in the file and a uniform probability over the "
                        "transitions. The value in the file gets weight "
                        "(1-loadsmoothtransition) while uniform gets weight "
                        "loadsmoothtransition. The default value is 0.5.")
        # [-x maxseconds]
        self.add_option('x', int, optional = True,
                        description = "This parameter specifies the maximum "
                        "number of seconds that can be spent optimizing the "
                        "model parameters. If it is less than 0, then there "
                        "is no limit and termination is based on maximum "
                        "number of iterations or a log likelihood change "
                        "criteria. The default value of this parameter is -1.")
        # [-z zerotransitionpower]
        self.add_option('z', int, optional = True,
                        description = "This parameter determines the threshold "
                        "at which to set extremely low transition "
                        "probabilities to 0 durining training. Setting "
                        "extremely low transition probabilities makes model "
                        "learning more efficient with essentially no impact "
                        "on the final results. If a transition probability "
                        "falls below 10^-zerotransitionpower during training "
                        "it is set to 0. Making this parameter to low and thus "
                        "the cutoff too high can potentially cause some "
                        "numerical instability. By default this parameter is "
                        "set to 8.")

    def runs(self, run_ids_connections_files):

        options = ['b', 'd', 'e', 'h', 'holdcolumnorder', 'init', 'm', 'nobed',
                   'nobrowser', 'noenrich', 'r', 's', 'stateordering', 't',
                   'x', 'z']

        set_options = [option for option in options if \
                       self.is_option_set_in_config(option)]

        option_list = list()
        for option in set_options:
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    option_list.append('-%s' % option)
                else:
                    option_list.append('-%s' % option)
            else:
                option_list.append('-%s' % option)
                option_list.append(str(self.get_option(option)))


        # TODO: 1. Extract the binary files into a directory
        # TODO: 2. 

        # We need to assemble the ChromHMM command that's all!
        for run_id in run_ids_connections_files.keys():
            input_paths = run_ids_connections_files[run_id]\
                          ['in/chromhmm_binarization']
            # Every file in input_paths should be in the same directory but
            # let's check:
            input_dir = input_paths[0]
            for f in input_paths[1:]:
                if not input_dir == os.path.dirname(f):
                    logger.error("Two binarization files from the same run ID"
                                 "are not in the same directory: %s and %s" %
                                 (input_paths[0], f))
                    sys.exit(1)

            with self.declare_run(run_id) as run:
                chromhmm = [self.get_tool('ChromHMM'),
                            'LearnModel']
                chromhmm.extend(option_list)
                chromhmm.append(input_dir)
                chromhmm.append(run.get_output_directory_du_jour_placeholder())
                chromhmm.append(self.get_option('numstates'))
                chromhmm.append(self.get_option('assembly'))
                
                run.add_output_file(
                    '',
                    ,
                    input_paths)
