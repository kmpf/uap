import os
uap_path = os.path.dirname(os.path.realpath(__file__))
activate_this_file = '%s/python_env/bin/activate_this.py' % uap_path
execfile(activate_this_file, dict(__file__=activate_this_file))
import argparse
import logging
import string
import sys
import yaml

'''
Adjust sys.path so everything we need can be found
'''

if uap_path not in sys.path:
    sys.path.append(uap_path)
include_path = '%s/include' % uap_path
if include_path not in sys.path:
    sys.path.append(include_path)
sources_path = '%s/include/sources' % uap_path
if sources_path not in sys.path:
    sys.path.append(sources_path)
steps_path = '%s/include/steps' % uap_path
if steps_path not in sys.path:
    sys.path.append(steps_path)
subcommand_path = '%s/include/subcommands' % uap_path
if subcommand_path not in sys.path:
    sys.path.append(subcommand_path)

from include.subcommands import *
from uaperrors import UAPError

def main():
    '''
    Set environment variables for git, so it does use the correct repository
    '''
    os.environ["GIT_DIR"] = "%s/.git" % uap_path
    os.environ["GIT_WORK_TREE"] = uap_path

    '''
    This script allows access to all commands which are provided by the pipeline.
    '''

    # Definition of common parser(s)

    common_parser = argparse.ArgumentParser(
        add_help=False,
        formatter_class=argparse.RawTextHelpFormatter,
        prog='uap'
    )

    common_parser.add_argument(
        "--even-if-dirty",
        dest="even_if_dirty",
        action="store_true",
        default=False,
        help="This option must be set if the local git repository "
        "contains uncommited changes.\n"
        "Otherwise uap will not run.")

    # Definition of the final parser

    parser = argparse.ArgumentParser(
        description="This script starts and controls analysis for 'uap'.",
        epilog="For complete documentation see: "
        "http://uap.readthedocs.org/en/latest/\n"
        "For citation use: ...\n"
        "For source code see: https://github.com/yigbt/uap",
        formatter_class=argparse.RawTextHelpFormatter,
        prog='uap'
        )

    parser.add_argument(
        "config",
        help="Path to YAML file that contains the pipeline configuration.\n"
        "The content of that file needs to follow the documentation.",
        metavar="<project-config>.yaml",
        nargs='?',
        type=argparse.FileType('r'))

    parser.add_argument(
        "-v", "--verbose",
        dest="verbose",
        action="count",
        default=1,
        help="Increase output verbosity")

    parser.add_argument(
        "--path",
        action="version",
        version=uap_path,
        help="Report the path of the UAP installation and exit.")

    parser.add_argument(
        "--debugging",
        dest="debugging",
        action="store_true",
        default=False,
        help="Print traceback on UAPError.")

    parser.add_argument(
        "--profiling",
        dest="profiling",
        action="store_true",
        default=False,
        help="Enable profiling save report in uap.cprof.")

    parser.add_argument(
        "--version",
        dest="version",
        action="version",
        version="%(prog)s 1.1",
        help = "Display version information.")

    subparsers = parser.add_subparsers(
        title="subcommands",
        description="Available subcommands.")

    '''
    The argument parser for 'fix-problems.py' is created here."
    '''

    fix_problems_parser = subparsers.add_parser(
        "fix-problems",
        help="Fixes problematic states by removing stall files.",
        description="",
        formatter_class=argparse.RawTextHelpFormatter,
        parents=[common_parser])

    fix_problems_parser.add_argument(
        "--cluster",
        dest="cluster",
        type=str,
        default="auto",
        help="Specify the cluster type. Default: [auto].")

    fix_problems_parser.add_argument(
        "--details",
        dest="details",
        action="store_true",
        default=False,
        help="Displays information about the files causing problems.")

    fix_problems_parser.add_argument(
        "--srsly",
        dest="srsly",
        action="store_true",
        default=False,
        help="Delete problematic files.")

    fix_problems_parser.set_defaults(func=fix_problems.main)

    '''
    The argument parser for 'render.py' is created here."
    '''

    render_parser = subparsers.add_parser(
        "render",
        help="Renders DOT-graphs displaying information of the analysis.",
        description="'render' generates DOT-graphs. Without arguments\n"
        "it takes the annotation file of each run and generates a graph,\n"
        "showing details of the computation.",
        formatter_class=argparse.RawTextHelpFormatter,
        parents=[common_parser])

    render_parser.add_argument(
        "--files",
        dest="files",
        action="store_true",
        default=False,
        help="Renders a graph showing all files of the analysis. "
        "[Not implemented yet!]")

    render_parser.add_argument(
        "--steps",
        dest="steps",
        action="store_true",
        default=False,
        help="Renders a graph showing all steps of the analysis and their "
        "connections.")

    render_parser.add_argument(
        "--simple",
        dest="simple",
        action="store_true",
        default=False,
        help="Simplify rendered graphs.")

    render_parser.add_argument(
        "--orientation",
        choices=['left-to-right', 'right-to-left', 'top-to-bottom'],
        dest="orientation",
        default='top-to-bottom',
        help="Defines orientation of the graph. Default: 'top-to-bottom'",
        type=str
    )

    render_parser.add_argument(
        "run",
        nargs='*',
        default=list(),
        type=str,
        help="Render only graphs for these runs.")

    render_parser.set_defaults(func=render.main)

#    '''
#    The argument parser for 'report.py' is created here."
#    '''
#
#    report_parser = subparsers.add_parser(
#        "report",
#        help="Generates reports of steps which can do so.",
#        description="",
#        formatter_class=argparse.RawTextHelpFormatter,
#        parents=[common_parser])
#
#    report_parser.add_argument(
#        "run",
#        nargs='*',
#        default=list(),
#        type=str,
#        help="Can take multiple step names as input. A step name is the name "
#        "of any entry in the 'steps:' section as defined in "
#        "'<project-config>.yaml'. To get a list of all run IDs please "
#        "run 'uap <project-config>.yaml status'.")
#
#    report_parser.set_defaults(func=report.main)

    '''
    The argument parser for 'run_locally.py' is created here."
    '''

    run_locally_parser = subparsers.add_parser(
        "run-locally",
        help="Executes the analysis on the local machine.",
        description="This command  starts 'uap' on the local machine. "
        "It can be used to start:\n"
        " * all runs of the pipeline as configured in <project-config>.yaml\n"
        " * all runs defined by a specific step in <project-config>.yaml\n"
        " * one or more steps\n"
        "To start the complete pipeline as configured in <project-config>.yaml "
        "execute:\n"
        "$ uap <project-config>.yaml run-locally\n"
        "To start a specific step execute:\n"
        "$ uap <project-config>.yaml run-locally <step_name>\n"
        "To start a specific run execute:\n"
        "$ uap <project-config>.yaml run-locally <step/run>\n"
        "The step_name is the name of an entry in the 'steps:' section "
        "as defined in '<project-config>.yaml'. A specific run is defined via "
        "its run ID 'step/run'. To get a list of all run IDs please run:\n"
        "$ uap <project-config>.yaml status",
        formatter_class=argparse.RawTextHelpFormatter,
        parents=[common_parser])

    run_locally_parser.add_argument(
        "run",
        nargs='*',
        default=list(),
        type=str,
        help="These runs are processed on the local machine.")

    run_locally_parser.set_defaults(func=run_locally.main)

    '''
    The argument parser for 'status.py' is created here.
    '''

    status_parser = subparsers.add_parser(
        "status",
        help="Displays information about the status of the analysis.",
        description="This script displays by default information about all "
        "runs of the pipeline as configured in '<project-config>.yaml'. But "
        "the displayed information can be narrowed down via command line "
        "options.\n"
        "IMPORTANT: Hints given by this script are just valid if the jobs were "
        "submitted to the cluster.",
        formatter_class=argparse.RawTextHelpFormatter,
        parents=[common_parser])

    status_parser.add_argument(
        "--cluster",
        dest="cluster",
        type=str,
        default="auto",
        help="Specify the cluster type. Default: [auto].")

    status_parser.add_argument(
        "--job-ids",
        dest="job_ids",
        action="store_true",
        default=False,
        help="Prints space seperated cluster job ids of all submitted jobs.")

    status_parser.add_argument(
        "--summarize",
        dest="summarize",
        action="store_true",
        default=False,
        help="Displays summarized information of the analysis.")

    status_parser.add_argument(
        "--graph",
        dest="graph",
        action="store_true",
        default=False,
        help="Displays the dependency graph of the analysis.")

    status_parser.add_argument(
        "--sources",
        dest="sources",
        action="store_true",
        default=False,
        help="Displays only information about the source runs.")

    status_parser.add_argument(
        "-r","--run",
        dest="run",
        nargs='*',
        default=list(),
        type=str,
        help="The status of these runs are displayed.")

    status_parser.set_defaults(func=status.main)

    '''
    The argument parser for 'steps.py' is created here.
    '''

    steps_parser = subparsers.add_parser(
        "steps",
        help="Displays information about the steps available in uap.",
        description="This script displays by default a list of all "
        "steps the pipeline can use.\n",
        formatter_class=argparse.RawTextHelpFormatter,
        parents=[common_parser])

    steps_parser.add_argument(
        "--details",
        dest="details",
        action="store_true",
        default=False,
        help="Displays dexription per step.")

    steps_parser.add_argument(
        "--show",
        dest="step",
        type=str,
        default="",
        help="Show the details of a specific step.")

    steps_parser.set_defaults(func=steps.main)

    '''
    The argument parser for 'submit-to-cluster.py' is created here."
    '''

    submit_to_cluster_parser = subparsers.add_parser(
        "submit-to-cluster",
        help="Submits the jobs created by uap to a cluster",
        description="This script submits all runs configured in "
        "<project-config>.yaml to a cluster. "
        "The configuration for the available cluster types is stored at "
        "/<path-to-uap>/cluster/cluster-specific-commands.yaml. "
        "The list of runs can be narrowed down to specific steps. "
        "All runs of the specified step will be submitted to the cluster. "
        "Also, individual runs IDs (step/run) can be used for submission.",
        formatter_class=argparse.RawTextHelpFormatter,
        parents=[common_parser])

    submit_to_cluster_parser.add_argument(
        "--cluster",
        dest="cluster",
        type=str,
        default="auto",
        help="Specify the cluster type. Default: [auto].")

    submit_to_cluster_parser.add_argument(
        "--first-error",
        dest="first_error",
        action="store_true",
        default=False,
        help="Print stderr of the first failed cluster job.")

    submit_to_cluster_parser.add_argument(
        "--job-ids",
        dest="job_ids",
        action="store_true",
        default=False,
        help="Prints space seperated cluster job ids of all submitted jobs.")

    submit_to_cluster_parser.add_argument(
        "--legacy",
        dest="legacy",
        action="store_true",
        default=False,
        help="Use none array cluster submission.")

    submit_to_cluster_parser.add_argument(
        "run",
        nargs='*',
        default=list(),
        type=str,
        help="Submit only these runs to the cluster.")

    submit_to_cluster_parser.set_defaults(func=submit_to_cluster.main)

    '''
    The argument parser for 'run-info.py' is created here."
    '''

    run_info_parser = subparsers.add_parser(
        "run-info",
        help="Displays information about certain source or processing runs.",
        description="",
        formatter_class=argparse.RawTextHelpFormatter,
        parents=[common_parser])

    run_info_parser.add_argument(
         "--sources",
        dest="sources",
        action="store_true",
        default=False,
        help="Displays only information about the source runs.")

    run_info_parser.add_argument(
        "-r","--run",
        dest="run",
        nargs='*',
        default=list(),
        type=str,
        help="Display run-info for these runs.")

    run_info_parser.set_defaults(func=run_info.main)

    '''
    The argument parser for 'volatilize.py' is created here."
    '''

    volatilize_parser = subparsers.add_parser(
        "volatilize",
        help="Saves disk space by volatilizing intermediate results",
        description="Save disk space by volatilizing intermediate results. "
        "Only steps marked with '_volatile: True' are considered.",
        formatter_class=argparse.RawTextHelpFormatter,
        parents=[common_parser])

    volatilize_parser.add_argument(
        "--details",
        dest="details",
        action="store_true",
        default=False,
        help="Shows which files can be volatilized.")

    volatilize_parser.add_argument(
        "--srsly",
        dest="srsly",
        action="store_true",
        default=False,
        help="Replaces files marked for volatilization with a placeholder.")

    volatilize_parser.set_defaults(func=volatilize.main)

    '''
    The argument parser for 'runtime_info.py' is created here."
    '''

    runtime_info_parser = subparsers.add_parser(
        "runtime-info",
        help="Provides Information about the runtime",
        description="Provides Information about the runtime.",
        formatter_class=argparse.RawTextHelpFormatter,
        parents=[common_parser])

    runtime_info_parser.add_argument(
        "--to-screen",
        dest="to_screen",
        action="store_true",
        default=False,
        help="Displays information to the screen instead of writing to file.")

    runtime_info_parser.add_argument(
        "--run_id",
        dest="run_id",
        action="store_true",
        default=False,
        help="Get information for a specific run-id.")

    runtime_info_parser.add_argument(
        "--file-path",
        dest="file_path",
        type=str,
        default='runtime_info.csv',
        help="Path of the output file. [runtime_info.csv]")

    runtime_info_parser.set_defaults(func=runtime_info.main)

    # get arguments and call the appropriate function
    args = parser.parse_args()
    # Add the path to this very file
    args.uap_path = uap_path

    # create logger object
    logger = _configure_logger(args.verbose)
    # call subcommand
    try:
        if args.profiling is True:
            import cProfile
            cProfile.runctx('args.func(args)', {'args': args}, {},
                    filename='uap.cprof')
        else:
            args.func(args)
    except UAPError as e:
        if args.debugging is True:
            raise
        else:
            sys.exit(1)

def _configure_logger(verbosity):
    logger = logging.getLogger("uap_logger")

    # create formatter for different log level
    debug_formatter = logging.Formatter(
        fmt = '[uap][%(levelname)s] %(funcName)s in %(filename)s: %(message)s '
    )
    info_formatter = logging.Formatter(
        fmt = '[uap][%(levelname)s]: %(message)s '
    )
    # create console handler
    ch = logging.StreamHandler()
    # set handler logging level
    ch.setLevel(logging.NOTSET)

    # Instantiate logger
    if verbosity == 0:
        # add formatter to ch
        ch.setFormatter(info_formatter)
        # set logger logging level
        logger.setLevel(logging.ERROR)
        logger.info("[uap] Set log level to ERROR\n")
    elif verbosity == 1:
        # add formatter to ch
        ch.setFormatter(info_formatter)
        # set logger logging level
        logger.setLevel(logging.WARNING)
        logger.info("[uap] Set log level to WARNING\n")
    elif verbosity == 2:
        # add formatter to ch
        ch.setFormatter(info_formatter)
        # set logger logging level
        logger.setLevel(logging.INFO)
        logger.info("[uap] Set log level to INFO\n")
    elif verbosity >= 3:
        # add formatter to ch
        ch.setFormatter(debug_formatter)
        # set logger logging level
        logger.setLevel(logging.DEBUG)
        logger.info("[uap] Set log level to DEBUG\n")

    sys.stderr.flush()

    # add ch to logger
    logger.addHandler(ch)

    return logger



if __name__ == '__main__':
    main()
