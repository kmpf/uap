..
  This is the documentation for uap. Please keep lines under 80 characters if
  you can and start each sentence on a new line as it decreases maintenance
  and makes diffs more readable.

.. title:: Command-Line Usage of uap

..
  This document aims to describe how to use **uap** via the command-line.

.. _cli_usage_uap:
#############################
Command-Line Usage of **uap**
#############################

**uap** uses Python's |argparse_link|.
Therefore, **uap** provides help infromation on the command-line::

  $ uap -h
  usage: uap [-h] [-v] [--version]
             [<project-config>.yaml]
             {fix-problems,render,report,run-locally,status,steps,
             submit-to-cluster,task-info,volatilize}
             ...

  This script starts and controls 'uap' analysis.

  positional arguments:
    <project-config>.yaml
                          Path to YAML file which holds the pipeline
                          configuration. It has to follow the structure given
                          in the documentation.

  optional arguments:
    -h, --help            show this help message and exit
    -v, --verbose         Increase output verbosity
    --version             Display version information.

  subcommands:
    Available subcommands.

    {fix-problems,render,report,run-locally,status,steps,submit-to-cluster,
    task-info,volatilize}
      fix-problems        Fixes problematic states by removing stall files.
      render              Renders DOT-graphs displaying information of the
                          analysis.
      report              Generates reports of steps which can do so.
      run-locally         Executes the analysis on the local machine.
      status              Displays information about the status of the analysis.
      steps               Displays information about the steps available in UAP.
      submit-to-cluster   Submits the jobs created by the seq-pipeline to a
                          cluster
      task-info           Displays information about certain source or
                          processing tasks.
      volatilize          Saves disk space by volatilizing intermediate results

  For further information please visit http://uap.readthedocs.org/en/latest/
  For citation use ...

Almost all subcommands require a YAML configuration file (see
:doc:`configuration`) **except** for ``uap steps``, which does not dependent
on a concrete analysis.

Everytime **uap** is started with a :ref:`configuration-of-uap` several things
happen:

1. The configuration file is read
2. The tools given in the :ref:`tools section <uap_config_tools_section>` are
   checked
3. The input files are checked
4. The state of all runs are calculated.

If any of these steps fails, **uap** will print an error message with and it
will crash.
This may seem a bit harsh, but after all, it's better to fail early than
to fail late if failing is unavoidable.

**uap** creates a symbolic link, if it does not exist already, pointing to the
:ref:`destination path <config-file-destination-path>` called
``<project-config>.yaml-out``.
The symbolic link is created in the directory containing the
``<project-config>.yaml``.

There are a couple of global command line parameters which are valid for all
scripts (well, actually, it's only one):

* ``--even-if-dirty``:
    Before doing anything else, the pipeline checks whether its source code
    has been modified in any way via Git.
    If yes, processing is stopped immediately unless this flag is specified.
    If you specify the flag, the fact that the repository was dirty will be
    recorded in all annotations which are produces *including* a full Git diff.

The subcommands are described in detail below.

.. _ExplanationOfSubcommands:
**************************
Explanation of Subcommands
**************************

.. _uap-steps:
steps
=====

The ``steps`` subcommand lists all available :ref:`source
<config_file_source_steps>` and :ref:`processing <config_file_processing_steps>`
steps::

  usage: uap [<project-config>.yaml] steps [-h] [--even-if-dirty] [--show STEP]

  This script displays by default a list of all steps the pipeline can use.

  optional arguments:
    -h, --help       show this help message and exit
    --even-if-dirty  Must be set if the local git repository contains uncommited
                     changes. Otherwise the pipeline will not start.
    --show STEP      Show the details of a specific step.

.. _uap-status:
status
======

The ``status`` subcommand lists all runs of an analysis.
A run is describes the concrete processing of a sample by a step.
Samples are usually defined at the source steps and are then propagated through
the analysis.
Here is the help message::

  $ uap <project-config>.yaml status -h
  usage: uap [<project-config>.yaml] status [-h] [--even-if-dirty]
                                            [--cluster CLUSTER] [--summarize]
                                            [--graph] [--sources]
                                            [-t [TASK [TASK ...]]]

  This script displays by default information about all tasks of the pipeline
  as configured in '<project-config>.yaml'. But the displayed information can
  be narrowed down via command line options.
  IMPORTANT: Hints given by this script are just valid if the jobs were
  submitted to the cluster.

  optional arguments:
    -h, --help            show this help message and exit
    --even-if-dirty       Must be set if the local git repository contains
                          uncommited changes. Otherwise the pipeline will not
                          start.
    --cluster CLUSTER     Specify the cluster type (sge, slurm), defaults to
                          auto.
    --summarize           Displays summarized information of the analysis.
    --graph               Displays the dependency graph of the analysis.
    --sources             Displays only information about the source runs.
    -t [TASK [TASK ...]], --task [TASK [TASK ...]]
                          Displays only the named task IDs. Can take multiple
                          task ID(s) as input. A task ID looks like this
                          'step_name/run_id'. A list of all task IDs is returned
                          by running:
                          $ uap <project-config>.yaml status


At any time, each run is in one of the following states:

* **waiting** -- the run is waiting for input files to appear, or its input
  files are not up-to-date regarding their respective dependencies
* **ready** -- all input files are present and up-to-date regarding their
  upstream input files (and so on, recursively), the run is ready and can
  be started
* **queued** -- the run is currently queued and will be started "soon"
  (only available if you use a compute cluster)
* **executing** -- the run is currently running on this or another machine
* **finished** -- all output files are in place and up-to-date



Here is an example output::

    $ uap <project-config>.yaml status
    Waiting tasks
    -------------
    [w] cufflinks/Sample_COPD_2023

    Ready tasks
    -----------
    [r] tophat2/Sample_COPD_2023

    Finished tasks
    --------------
    [f] cutadapt/Sample_COPD_2023-R1
    [f] cutadapt/Sample_COPD_2023-R2
    [f] fix_cutadapt/Sample_COPD_2023

    tasks: 5 total, 1 waiting, 1 ready, 3 finished

To get a more concise summary, specify ``--summarize``::

    $ uap <project-config>.yaml status --summarize
    Waiting tasks
    -------------
    [w]   1 cufflinks

    Ready tasks
    -----------
    [r]   1 tophat2

    Finished tasks
    --------------
    [f]   2 cutadapt
    [f]   1 fix_cutadapt

    tasks: 5 total, 1 waiting, 1 ready, 3 finished

...or print a fancy ASCII art graph with ``--graph``::

    $ uap <project-config>.yaml status --graph
    samples (1 finished)
    └─cutadapt (2 finished)
      └─fix_cutadapt (1 finished)
        └─tophat2 (1 ready)
          └─cufflinks (1 waiting)


Detailed information about a specific task can be obtained by specifying the
run ID on the command line::

  $ uap index_mycoplasma_genitalium_ASM2732v1_genome.yaml status -t \
    bowtie2_index/Mycoplasma_genitalium_index-download --even
  [uap] Set log level to ERROR
  output_directory: genomes/bacteria/Mycoplasma_genitalium/bowtie2_index/
                    Mycoplasma_genitalium_index-download-cMQPtBxs
  output_files:
    out/bowtie_index:
      Mycoplasma_genitalium_index-download.1.bt2: &id001
      - genomes/bacteria/Mycoplasma_genitalium/Mycoplasma_genitalium.ASM2732v1.fa
      Mycoplasma_genitalium_index-download.2.bt2: *id001
      Mycoplasma_genitalium_index-download.3.bt2: *id001
      Mycoplasma_genitalium_index-download.4.bt2: *id001
      Mycoplasma_genitalium_index-download.rev.1.bt2: *id001
      Mycoplasma_genitalium_index-download.rev.2.bt2: *id001
  private_info: {}
  public_info: {}
  run_id: Mycoplasma_genitalium_index-download
  state: FINISHED

This data structure is called the "run info" of a certain run and it
represents a kind of plan which includes information about which output
files will be generated and which input files they depend on -- this is
stored in ``output_files``.

Source steps can be viewed separately by specifying ``--sources``::

    $ uap <project-config>.yaml status --sources
    [uap] Set log level to ERROR
    M_genitalium_genome/download

.. _uap-task-info:
task-info
=========

The ``task-info`` subcommand writes the commands which were or will be executed
to the terminal in the form of a semi-functional BASH script.
Semi-functional means that at the moment output redirections for some commands
are missing in the BASH script.
Also included are the ``run info`` informations as already described for the
``status`` subcommand.

An example output showing the download of the *Mycoplasma genitalium* genome::

  $ uap index_mycoplasma_genitalium_ASM2732v1_genome.yaml task-info -t \
    M_genitalium_genome/download --even

    [uap] Set log level to ERROR
    #!/usr/bin/env bash

    # M_genitalium_genome/download -- Report
    # ======================================
    #
    # output_directory: genomes/bacteria/Mycoplasma_genitalium/M_genitalium_genome/download-5dych7Xj
    # output_files:
    #   out/raw:
    #     genomes/bacteria/Mycoplasma_genitalium/Mycoplasma_genitalium.ASM2732v1.fa: []
    # private_info: {}
    # public_info: {}
    # run_id: download
    # state: FINISHED
    #
    # M_genitalium_genome/download -- Commands
    # ========================================

    # 1. Group of Commands -- 1. Command
    # ----------------------------------

    curl ftp://ftp.ncbi.nih.gov/genomes/genbank/bacteria/Mycoplasma_genitalium/latest_assembly_versions/GCA_000027325.1_ASM2732v1/GCA_000027325.1_ASM2732v1_genomic.fna.gz

    # 2. Group of Commands -- 1. Command
    # ----------------------------------

    ../tools/compare_secure_hashes.py --algorithm md5 --secure-hash a3e6e5655e4996dc2d49f876be9d1c27 genomes/bacteria/Mycoplasma_genitalium/M_genitalium_genome/download-5dych7Xj/L9PXBmbPKlemghJGNM97JwVuzMdGCA_000027325.1_ASM2732v1_genomic.fna.gz

    # 3. Group of Commands -- 1. Pipeline
    # -----------------------------------

    pigz --decompress --stdout --processes 1 genomes/bacteria/Mycoplasma_genitalium/M_genitalium_genome/download-5dych7Xj/L9PXBmbPKlemghJGNM97JwVuzMdGCA_000027325.1_ASM2732v1_genomic.fna.gz | dd bs=4M of=/home/hubert/develop/uap/example-configurations/genomes/bacteria/Mycoplasma_genitalium/Mycoplasma_genitalium.ASM2732v1.fa

This subcommand allows the user to run parts of the analysis and manually control
for causes of failure.


.. _uap-run-locally:
run-locally
===========

The ``run-locally`` subcommand runs all non-finished runs (or a subset)
sequentially on the local machine.
Feel free to cancel this script at any time, it won't put your project in a
unstable state.
However, if the ``run-locally`` subcommand receives a SIGKILL signal, the
currently executing job will continue to run and the corresponding run
will be reported as ``executing`` by calling ``status`` subcommand for five more
minutes (SIGTERM should be fine and exit gracefully but *doesn't just yet*).
After that time, you will be warned that a job is marked as being currently
run but no activity has been seen for a while, along with further
instructions about what to do in such a case (don't worry, it shouldn't
happen by accident).

To execute one or more certain runs, specify the run IDs on the command
line.
To execute all runs of a certain step, specify the step name on the command
line.

This subcommands usage information::

  $ uap <project-config>.yaml run-locally -h
  usage: uap [<project-config>.yaml] run-locally [-h] [--even-if-dirty]
                                               [step_task [step_task ...]]

  This command  starts 'uap' on the local machine. It can be used to start:
  * all tasks of the pipeline as configured in <project-config>.yaml
  * all tasks defined by a specific step in <project-config>.yaml
  * one or more steps
  To start the complete pipeline as configured in <project-config>.yaml execute:
    $ uap <project-config>.yaml run-locally
  To start a specific step execute:
    $ uap <project-config>.yaml run-locally <step_name>
  To start a specific task execute:
    $ uap <project-config>.yaml run-locally <step_name/run_id>
  The step_name is the name of an entry in the 'steps:' section as defined in
  '<project-config>.yaml'. A specific task is defined via its task ID
  'step_name/run_id'. A list of all task IDs is returned by running:
    $ uap <project-config>.yaml status

  positional arguments:
    step_task        Can take multiple step names as input. A step name is the
                     name of any entry in the 'steps:' section as defined in
                     '<config>.yaml'. A list of all task IDs is returned by running:
                       $ uap <project-config>.yaml status.

  optional arguments:
    -h, --help       show this help message and exit
    --even-if-dirty  Must be set if the local git repository contains uncommited
                     changes. Otherwise the pipeline will not start.

.. NOTE:: Why is it safe to cancel the pipeline?
    The pipeline is written in a way which expects processes to fail or
    cluster jobs to disappear without notice.
    This problem is mitigated by a design which relies on file presence and
    file timestamps to determine whether a run is finished or not.
    Output files are automatically written to temporary locations and later
    moved to their real target directory, and it is not until the last file
    rename operation has finished that a run is regarded as finished.

.. _uap-submit-to-cluster:
submit-to-cluster
=================

The ``submit-to-cluster`` subcommand determines which runs still have to be
carried out and which supported cluster engine is available.
It then submits the jobs to the cluster if a cluster engine has been found.
Dependencies are passed to cluster engine in a way that jobs that depend on
other jobs won't get scheduled until their dependencies have been satisfied.
The files ``qsub-template.sh`` and ``sbatch-template.sh`` are used to submit
jobs.
Fields with ``#{ }`` are substituted with appropriate values.
Each submitted job calls **uap** with the ``run-locally`` subcommand on the
cluster nodes.

The file ``quotas.yaml`` can be used to define different quotas for different
systems:

.. code-block:: yaml

    "frontend[12]":
        default: 5
        cutadapt: 100

In the example above, a default quota of 5 is defined for hosts with a
hostname of ``frontend1`` or ``frontend2`` (the name is a regular expression).
A quota of 5 means that no more than 5 jobs of one kind will be run in
parallel.
Different quotas can be defined for each step: because ``cutadapt`` is
highly I/O-efficient, it has a higher quota.

This subcommand provides usage information::

  $ uap <project-config>.yaml submit-to-cluster -h
  usage: uap [<project-config>.yaml] submit-to-cluster [-h] [--even-if-dirty]
                                                       [--cluster CLUSTER]
                                                       [step_task [step_task ...]]

  This script submits all tasks configured in <project-config>.yaml to a
  SGE/OGE/UGE or SLURM cluster. The list of tasks can be narrowed down by
  specifying a step name (in which case all runs of this steps will be considered)
  or individual tasks (step_name/run_id).

  positional arguments:
    step_task          Can take multiple step names as input. A step name is
                       the name of any entry in the 'steps:' section as defined
                       in '<project-config>.yaml'. A list of all task IDs is
                       returned by running:
                         $ uap <project-config>.yaml status

    optional arguments:
      -h, --help         show this help message and exit
      --even-if-dirty    Must be set if the local git repository contains
                         uncommited changes. Otherwise the pipeline will not
                         start.
      --cluster CLUSTER  Specify the cluster type. Choices: [auto, sge, slurm].
                         Default: [auto].


.. _uap-fix-problems:
fix-problems
============

The ``fix-problems`` subcommand removes temporary files written by **uap** if
they are not required anymore.
This subcommand provides usage information::

  $ uap <project-config>.yaml fix-problems -h
  usage: uap [<project-config>.yaml] fix-problems [-h] [--even-if-dirty]
                                                  [--cluster CLUSTER]
                                                  [--details] [--srsly]

  optional arguments:
    -h, --help         show this help message and exit
    --even-if-dirty    Must be set if the local git repository contains
                       uncommited changes. Otherwise the pipeline will not start.
    --cluster CLUSTER  Specify the cluster type (sge, slurm), defaults to auto.
    --details          Displays information about problematic files which need
                       to be deleted to fix problem.
    --srsly            Deletes problematic files.


**uap** writes temporary files to indicate if a job is queued or executed.
Sometimes (especially on the compute cluster) jobs fail, without even starting
**uap**.
This leaves the temporary file, written on job submission, indicating that a run
was queued on the cluster without process (because it already failed).
The ``status`` subcommand will inform the user if ``fix-problems`` needs to be
executed to clean up the mess.
The hint given by ``status`` would look like::

  Warning: There are 10 tasks marked as queued, but they do not seem to be queued
  Hint: Run 'uap <project-config>.yaml fix-problems --details' to see the details.
  Hint: Run 'uap <project-config>.yaml fix-problems --srsly' to fix these problems
        (that is, delete all problematic ping files).

Be nice and do as you've told.
Now you are able to resubmit your runs to the cluster.
You've fixed the problem, haven't you?

.. _uap-volatilize:
volatilize
==========

The ``volatilize`` subcommand is useful to reduce the required disk space of
your analysis.
It works only in conjunction with the :ref:`_volatile <uap-volatile>` keyword
set in the :ref:`configuration file <configuration_of_uap>`.
As already mentioned there, steps marked as ``_volatile`` compute their output
files as normal but they can be deleted if their dependent steps are finished.

This subcommand provides usage information::

  $ uap <project-config>.yaml volatilize -h
  usage: uap [<project-config>.yaml] volatilize [-h] [--even-if-dirty]
                                                [--details] [--srsly]

  Save disk space by volatilizing intermediate results. Only steps marked with
  '_volatile: True' are considered.

  optional arguments:
    -h, --help       show this help message and exit
    --even-if-dirty  Must be set if the local git repository contains uncommited
                     changes. Otherwise the pipeline will not start.
    --details        Shows which files can be volatilized.
    --srsly          Replaces files marked for volatilization with a placeholder.

After running ``volatilize --srsly`` the output files of the volatilized step
are replaced by placeholder files.
The placeholder files have the same name as the original files suffixed with
``.volatile.placeholder.yaml``.

.. _uap-render:
render
======

The ``render`` subcommand generates graphs using graphviz.
The graphs either show the complete analysis or the execution of a single run.
At the moment ``--simple`` only has an effect in combination with ``--steps``.

This subcommand provides usage information::

   $ uap <project-config>.yaml render -h

   usage: uap [<project-config>.yaml] render [-h] [--even-if-dirty] [--files]
                                             [--steps] [--simple]
                                             [step_task [step_task ...]]

   'render' generates DOT-graphs. Without arguments it takes the log file of
   each task and generates a graph, showing details of the computation.

   positional arguments:
     step_task        Displays only the named task IDs. Can take multiple task
                      ID(s) as input. A task ID looks like this
                      'step_name/run_id'. A list of all task IDs is returned by
                      running 'uap <project-config>.yaml status'.

   optional arguments:
     -h, --help       show this help message and exit
     --even-if-dirty  Must be set if the local git repository contains
                      uncommited changes. Otherwise the pipeline will not start.
     --files          Renders a graph showing all files of the analysis.
                      [Not implemented yet!]
     --steps          Renders a graph showing all steps of the analysis and their
                      connections.
     --simple         Rendered graphs are simplified.



.. _uap-report:
report
======

The ``report`` subcommand is at the moment experimental.
It might be used to create standardized output to enable easy loading and
processing in downstream tools e.g. ``R``.

.. |argparse_link| raw:: html

   <a href="https://docs.python.org/2.7/library/argparse.html" target="_blank">argparse</a>.
