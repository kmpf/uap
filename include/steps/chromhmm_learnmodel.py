import os
from logging import getLogger
import tarfile

from abstract_step import AbstractStep

logger = getLogger('uap_logger')

class ChromHmmLearnModel(AbstractStep):
    '''
    This command takes a directory with a set of binarized data files and learns
    a chromatin state model. Binarized data files have "_binary" in the file
    name. The format for the binarized data files are that the first line
    contains the name of the cell separated by a tab with the name of the
    chromosome. The second line contains in tab delimited form the name of each
    mark. The remaining lines correspond to consecutive bins on the chromosome.
    The remaining lines in tab delimited form corresponding to each mark, with a
    "1" for a present call or "0" for an absent call and a "2" if the data is
    considered missing at that interval for the mark.
    '''

    def __init__(self, pipeline):
        super(ChromHmmLearnModel, self).__init__(pipeline)

        self.set_cores(8)
        
        self.add_connection('in/cellmarkfiletable')
        self.add_connection('in/chromhmm_binarization')
        self.add_connection('out/chromhmm_model')
        
        self.require_tool('ChromHMM')
        self.require_tool('ls')
        self.require_tool('mkdir')
        self.require_tool('rm')
        self.require_tool('tar')
        self.require_tool('xargs')
        
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
                        "2012). 'random' - randomly initializes the parameters "
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

        options = ['b', 'color', 'd', 'e', 'h', 'holdcolumnorder', 'init', 'l',
                   'm', 'nobed', 'nobrowser', 'noenrich',
                   # 'printposterior', 'printstatebyline',
                   'r', 's', 'stateordering', 't', 'x', 'z']
        file_options = ['assembly', 'l', 'm']

        set_options = [option for option in options if \
                       self.is_option_set_in_config(option)]

        option_list = list()
        for option in set_options:
            if isinstance(self.get_option(option), bool):
                # Only set option if it is True
                if self.get_option(option):
                    option_list.append('-%s' % option)
            else:
                value = str(self.get_option(option))
                if option in file_options:
                    value = os.path.abspath(value)
                option_list.append('-%s' % option)
                option_list.append(value)

        for run_id in run_ids_connections_files.keys():
            # The input_paths should be a single tar.gz file
            input_paths = run_ids_connections_files[run_id]\
                          ['in/chromhmm_binarization']
            # Test the input_paths (at least a bit)
            if len(input_paths) != 1 or not input_paths[0].endswith('.tar.gz'):
                logger.error("Expected single tar.gz file via "
                             "'in/chromhmm_binarization' for run %s, but got "
                             "this %s" % (run_id, ", ".join(input_paths)))
                sys.exit(1)


            # read tar file and get names of included files
#            with tarfile.open(name = input_paths[0], mode = 'r:gz') as tar:
#                tar.list()

            with self.declare_run(run_id) as run:
                with run.new_exec_group() as pre_chromhmm:
                    # 1. Extract the binary files into a directory
                    # 1.1 Get name of temporary input directory
                    input_dir =  run.add_temporary_directory(
                        '%s_binary_files' % run_id)
                    # 1.2 Create temporary input directory
                    mkdir = [self.get_tool('mkdir'), input_dir]
                    pre_chromhmm.add_command(mkdir)
                    # 1.3 Extract the binary files into temporary input directory
                    tar = [self.get_tool('tar'),
                           '--extract',
                           '--gzip',
                           '--verbose',
                           '--directory',
                           input_dir,
                           '--file',
                           input_paths[0] ]
                    pre_chromhmm.add_command(tar)
                    # 1.4 Get name of temporary output directory
                    output_dir = run.add_temporary_directory(
                        '%s_chromhmm_model' % run_id)
                    # 1.5 Create temporary output directory
                    mkdir = [self.get_tool('mkdir'), output_dir]
                    pre_chromhmm.add_command(mkdir)

                with run.new_exec_group() as learnmodel:
                    # 2. Assemble ChromHMM LearnModel command
                    chromhmm = [self.get_tool('ChromHMM'),
                                'LearnModel']
                    chromhmm.extend(option_list)
                    chromhmm.append(input_dir)
                    chromhmm.append(output_dir)
                    chromhmm.append(str(self.get_option('numstates')))
                    chromhmm.append(os.path.abspath(self.get_option('assembly')))
                    learnmodel.add_command(chromhmm)

                with run.new_exec_group() as pack_model:
                    # 3. Pack the output files of ChromHMM LearnModel
                    with pack_model.add_pipeline() as pack_model_pipe:
                        # 3.1 List content of output directory
                        ls = [self.get_tool('ls'), '-1', output_dir]
                        # 3.2 Pipe ls output
                        pack_model_pipe.add_command(ls)
                        # 3.3 Use xargs to call tar (circumventing glob pattern)
                        xargs = [self.get_tool('xargs'),
                                 '--delimiter', '\n',
                                 self.get_tool('tar'),
                                 '--create',
                                 '--directory',
                                 output_dir,
                                 '--gzip',
                                 '--remove-files',
                                 '--verbose',
                                 '--file',
                                 run.add_output_file(
                                     'chromhmm_model',
                                     '%s_model_files.tar.gz' % run_id,
                                     input_paths)]
                        pack_model_pipe.add_command(xargs)

                with run.new_exec_group() as rm_binary_files:
                    # 4. Remove the unpacked binary files
                    with rm_binary_files.add_pipeline() as rm_binary_pipe:
                        # 4.1 List content of output directory
                        ls = [self.get_tool('ls'), '-1', input_dir]
                        # 4.2 Pipe ls output
                        rm_binary_pipe.add_command(ls)
                        # 4.3 Use xargs to call tar (circumventing glob pattern)
                        xargs = [self.get_tool('xargs'),
                                 '--delimiter', '\n',
                                 '-I', '*',
                                 self.get_tool('rm'),
                                 '--verbose',
                                 os.path.join(input_dir, '*')
                        ]
                        rm_binary_pipe.add_command(xargs)
