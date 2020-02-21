..
  This is the documentation for uap. Please keep lines under 80 characters if
  you can and start each sentence on a new line as it decreases maintenance
  and makes diffs more readable.

.. title:: Command-Line Usage of uap

..
  This document aims to describe how to use **uap** via the command-line.

.. _cli_usage_uap:

*****************************
Command-Line Usage of **uap**
*****************************

**uap** uses Python's |argparse_link|.
Therefore, **uap** provides help information on the command-line::

    $ uap -h
    usage: uap [-h] [-v] [--path] [--debugging] [--profiling] [--version]
               [<project-config>.yaml]
               {fix-problems,render,run-locally,status,steps,submit-to-cluster,run-info,volatilize,runtime-info}
               ...

    This script starts and controls analysis for 'uap'.

    positional arguments:
      <project-config>.yaml
                            Path to YAML file that contains the pipeline configuration.
                            The content of that file needs to follow the documentation.

    optional arguments:
      -h, --help            show this help message and exit
      -v, --verbose         Increase output verbosity
      --path                Report the path of the UAP installation and exit.
      --debugging           Print traceback on UAPError.
      --profiling           Enable profiling save report in uap.cprof.
      --version             Display version information.

    subcommands:
      Available subcommands.

      {fix-problems,render,run-locally,status,steps,submit-to-cluster,run-info,volatilize,runtime-info}
        fix-problems        Fixes problematic states by removing stall files.
        render              Renders DOT-graphs displaying information of the analysis.
        run-locally         Executes the analysis on the local machine.
        status              Displays information about the status of the analysis.
        steps               Displays information about the steps available in uap.
        submit-to-cluster   Submits the jobs created by uap to a cluster
        run-info            Displays information about certain source or processing runs.
        volatilize          Saves disk space by volatilizing intermediate results
        runtime-info        Provides Information about the runtime

    For complete documentation see: http://uap.readthedocs.org/en/latest/
    For citation use: ...
    For source code see: https://github.com/yigbt/uap

Almost all subcommands require a YAML configuration file (see
:ref:`analysis_configuration`) **except** for ``uap steps``, which works
independent of an analysis configuration file.

Everytime **uap** is started with a
:ref:`analysis configuration file <analysis_configuration>` the following actions
happen:

1. Configuration file is read
2. Tools given in the :ref:`tools section <uap_config_tools_section>` are
   checked
3. Input files are checked
4. State of all runs are calculated

If any of these steps fail, **uap** will exit and print an error message.

**uap** will create a symbolic link, if it does not exist already, pointing to
the :ref:`destination path <config-file-destination-path>` called
``<project-config>.yaml-out``.
The symbolic link is created in the directory containing the
``<project-config>.yaml``.

There are a couple of global command line parameters which are valid for all
scripts:

``--even-if-dirty`` or short ``--even``:
    If this parameter appears **uap** will work even if uncommited changes
    to its source code are detected.
    **uap** would otherwise immediately stop working.
    If you specify this flag, the repositories state is recorded in all
    annotation files created by this process.
    A full Git diff is *included* as well.

``--no-tool-check``
    This disables the tool version check and can speed up, e.g., job
    submisstion with ``submit-to-cluster`` but is not recommended to
    use in ``status`` sice missing tool version may render a task
    state ``changed``.

``--debugging``
    If this parameter apperas **uap** will print a traceback on any
    error that occures instead of just the error message.

``--profiling``
    With this parameter the **uap** will use attempt to use cProfile
    and save any profiling in a file ``uap.cprof`` of the current
    user working directory.


.. _subcommands:

Subcommands
===========

Here an overview of all the available subcommands are given.

.. _uap-steps:

``steps`` Subcommand
--------------------

The ``steps`` subcommand lists all available :ref:`source
<config_file_source_steps>` and :ref:`processing <config_file_processing_steps>`
steps::

  $ uap steps -h
  usage: uap [<project-config>.yaml] steps [-h] [--even-if-dirty] [--show STEP]

  This script displays by default a list of all steps the pipeline can use.

  optional arguments:
    -h, --help       show this help message and exit
    --even-if-dirty  This option must be set if the local git repository
                     contains uncommited changes.
                     Otherwise uap will not run.
    --show STEP      Show the details of a specific step.


.. _uap-status:

``status`` Subcommand
---------------------

The ``status`` subcommand lists all runs of an analysis.
A run is describes the concrete processing of a sample by a step.
Samples are usually defined at the source steps and are then propagated through
the analysis.
Here is the help message::

  $ uap <project-config>.yaml status -h
  usage: uap [<project-config>.yaml] status [-h] [--even-if-dirty]
                                            [--no-tool-checks]
                                            [--cluster CLUSTER] [--details]
                                            [--job-ids] [--summarize] [--graph]
                                            [--hash] [--sources]
                                            [-r [RUN [RUN ...]]]

  This script displays by default information about all runs of the pipeline as configured in '<project-config>.yaml'. But the displayed information can be narrowed down via command line options.
  IMPORTANT: Hints given by this script are just valid if the jobs were submitted to the cluster.

  optional arguments:
    -h, --help            show this help message and exit
    --even-if-dirty       This option must be set if the local git repository contains uncommited changes.
                          Otherwise uap will not run.
    --no-tool-checks      This option disables the otherwise mandatory checks for tool availability and version
    --cluster CLUSTER     Specify the cluster type. Default: [auto].
    --details             Displays information about changed tasks.
    --job-ids             Prints space seperated cluster job ids of all submitted jobs.
    --summarize           Displays summarized information of the analysis.
    --graph               Displays the dependency graph of the analysis.
    --hash                Compare sha256sums of existing files with the logged values.
    --sources             Displays only information about the source runs.
    -r [RUN [RUN ...]], --run [RUN [RUN ...]]
                          The status of these runs are displayed.

At any time, each run is in one of the following states:

* ``[w]aiting`` -- the run is waiting for input files to appear, or its input
  files are not up-to-date regarding their respective dependencies
* ``[r]eady`` -- all input files are present and up-to-date regarding their
  upstream input files (and so on, recursively), the run is ready and can
  be started
* ``[q]ueued`` -- the run is currently queued and will be started "soon"
  (only available if you use a compute cluster)
* ``[e]xecuting`` -- the run is currently running on this or another machine
* ``[f]inished`` -- all output files are in place and up-to-date
* ``[c]hanged`` -- all output files are in place but the configuration,
  parent or the commands to execute changed
* ``[b]ad`` -- an error was caught during execution



Here is an example output::

    $ uap <project-config>.yaml status
    Waiting tasks
    -------------
    [w] fasta_index/download
    [w] segemehl_index/Mycoplasma_genitalium_genome-download

    Ready tasks
    -----------
    [r] bowtie2_index/Mycoplasma_genitalium_index-download
    [r] bwa_index/Mycoplasma_genitalium_index-download

    Finished tasks
    --------------
    [f] M_genitalium_genome/download

    tasks: 5 total, 2 waiting, 2 ready, 1 finished

To get a more concise summary, specify ``--summarize``::

    $ uap <project-config>.yaml status --summarize
    Waiting tasks
    -------------
    [w]   1 fasta_index
    [w]   1 segemehl_index

    Ready tasks
    -----------
    [r]   1 bowtie2_index
    [r]   1 bwa_index

    Finished tasks
    --------------
    [f]   1 M_genitalium_genome

    tasks: 5 total, 2 waiting, 2 ready, 1 finished

... or print a fancy ASCII art graph with ``--graph``::

    $ uap <project-config>.yaml status --graph
    M_genitalium_genome (raw_url_source) [1 finished]
    └─│─│─│─bowtie2_index (bowtie2_generate_index) [1 ready]
      └─│─│─bwa_index (bwa_generate_index) [1 ready]
        └─│─fasta_index (samtools_faidx) [1 waiting]
          └─segemehl_index (segemehl_generate_index) [1 waiting]

Detailed information about a specific task can be obtained by specifying the
run ID on the command line::

  $ uap index_mycoplasma_genitalium_ASM2732v1_genome.yaml status -r \
    bowtie2_index/Mycoplasma_genitalium_index-download
  output_directory: genomes/bacteria/Mycoplasma_genitalium/bowtie2_index/Mycoplasma_genitalium_index-download-ZsvbSjtK
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

This is the known data for run
``bowtie2_index/Mycoplasma_genitalium_index-download``.
It contains information about the output folder, the output files and the
input files they depend on as well as the run ID and the run state.

Source steps can be viewed separately by specifying ``--sources``::

    $ uap <project-config>.yaml status --sources
    M_genitalium_genome/download

.. _uap-run-info:

``run-info`` Subcommand
-----------------------

The ``run-info`` subcommand displays the commands issued for a given run.
The output looks like a BASH script, but might not be functional.
This is due to the fact that output redirections for some commands
are missing in the BASH script.
The output includes also the information as shown by the ``status -r <run-ID>``
subcommand.

An example output showing the download of the *Mycoplasma genitalium* genome::

  $ uap index_mycoplasma_genitalium_ASM2732v1_genome.yaml run-info --even -r M_genitalium_genome/download
  #!/usr/bin/env bash

  # M_genitalium_genome/download -- Report
  # ======================================
  #
  # output_directory: genomes/bacteria/Mycoplasma_genitalium/M_genitalium_genome/download-7RncJ4tr
  # output_files:
  #   out/raw:
  #     genomes/bacteria/Mycoplasma_genitalium/Mycoplasma_genitalium.ASM2732v1.fa: []
  # private_info: {}
  # public_info: {}
  # run_id: download
  #
  # M_genitalium_genome/download -- Commands
  # ========================================

  # 1. Group of Commands -- 1. Command
  # ----------------------------------

  curl ftp://ftp.ncbi.nih.gov/genomes/genbank/bacteria/Mycoplasma_genitalium/latest_assembly_versions/GCA_000027325.1_ASM2732v1/GCA_000027325.1_ASM2732v1_genomic.fna.gz

  # 2. Group of Commands -- 1. Command
  # ----------------------------------

  ../tools/compare_secure_hashes.py --algorithm md5 --secure-hash f02c78b5f9e756031eeaa51531517f24 genomes/bacteria/Mycoplasma_genitalium/M_genitalium_genome/download-7RncJ4tr/L9PXBmbPKlemghJGNM97JwVuzMdGCA_000027325.1_ASM2732v1_genomic.fna.gz

  # 3. Group of Commands -- 1. Pipeline
  # -----------------------------------

  pigz --decompress --stdout --processes 1 genomes/bacteria/Mycoplasma_genitalium/M_genitalium_genome/download-7RncJ4tr/L9PXBmbPKlemghJGNM97JwVuzMdGCA_000027325.1_ASM2732v1_genomic.fna.gz | dd bs=4M of=/home/hubert/develop/uap/example-configurations/genomes/bacteria/Mycoplasma_genitalium/Mycoplasma_genitalium.ASM2732v1.fa


This subcommand enables the user to manually run parts of the analysis without
**uap**.
That can be helpful for debugging steps during development.

.. _uap-run-locally:

``run-locally`` Subcommand
--------------------------

The ``run-locally`` subcommand runs all non-finished runs (or a specified
subset) sequentially on the local machine.
The execution can be cancelled at any time, it won't put your project in a
unstable state.
However, if the ``run-locally`` subcommand receives a |sigkill_link| signal, the
currently executing job will continue to run and the corresponding run
will be reported as ``executing`` by calling ``status`` subcommand for five more
minutes (|sigterm_link| should be fine and exit gracefully but
*doesn't just yet*).
After that time, you will be warned that a job is marked as being currently
run but no activity has been seen for a while, along with further
instructions about what to do in such a case (don't worry, it shouldn't
happen by accident).

Specify a set of run IDs to execute only those runs.
Specify the name of a step to execute all ready runs of that step.

This subcommands usage information::

  $ uap index_mycoplasma_genitalium_ASM2732v1_genome.yaml run-locally -h
  usage: uap [<project-config>.yaml] run-locally [-h] [--even-if-dirty]
                                                 [--no-tool-checks] [--force]
                                                 [--ignore]
                                                 [run [run ...]]

  This command  starts 'uap' on the local machine. It can be used to start:
   * all runs of the pipeline as configured in <project-config>.yaml
   * all runs defined by a specific step in <project-config>.yaml
   * one or more steps
  To start the complete pipeline as configured in <project-config>.yaml execute:
  $ uap <project-config>.yaml run-locally
  To start a specific step execute:
  $ uap <project-config>.yaml run-locally <step_name>
  To start a specific run execute:
  $ uap <project-config>.yaml run-locally <step/run>
  The step_name is the name of an entry in the 'steps:' section as defined in '<project-config>.yaml'. A specific run is defined via its run ID 'step/run'. To get a list of all run IDs please run:
  $ uap <project-config>.yaml status

  positional arguments:
    run               These runs are processed on the local machine.

  optional arguments:
    -h, --help        show this help message and exit
    --even-if-dirty   This option must be set if the local git repository contains uncommited changes.
                      Otherwise uap will not run.
    --no-tool-checks  This option disables the otherwise mandatory checks for tool availability and version
    --force           Force to overwrite changed tasks.
    --ignore          Ignore chages of tasks and consider them finished.

.. NOTE:: Why is it safe to cancel the pipeline?
    The pipeline is written in a way which expects processes to fail or
    cluster jobs to disappear without notice.
    This problem is mitigated by a design which relies on file presence and
    file timestamps to determine whether a run is finished or not.
    Output files are automatically written to temporary locations and later
    moved to their real target directory, and it is not until the last file
    rename operation has finished that a run is regarded as finished.

.. _uap-submit-to-cluster:

``submit-to-cluster`` Subcommand
--------------------------------

The ``submit-to-cluster`` subcommand determines which runs still need to be
executed and which supported cluster engine is available.
It submits a job for every run to the cluster if a cluster engine could be
detected.
Dependencies are passed to cluster engine in a way that jobs that depend on
other jobs won't get scheduled until their dependencies have been satisfied.
For more information read about the
:ref:`cluster configuration <cluster_configuration>` and the
:ref:`submit script template <submit_template>`.
Each submitted job calls **uap** with the ``run-locally`` subcommand on the
executing cluster node.

Here is the usage information::

  $ uap index_mycoplasma_genitalium_ASM2732v1_genome.yaml submit-to-cluster -h
  usage: uap [<project-config>.yaml] submit-to-cluster [-h] [--even-if-dirty]
                                                       [--no-tool-checks]
                                                       [--cluster CLUSTER]
                                                       [--legacy] [--force]
                                                       [--ignore]
                                                       [run [run ...]]

  This script submits all runs configured in <project-config>.yaml to a cluster. The configuration for the available cluster types is stored at /<path-to-uap>/cluster/cluster-specific-commands.yaml. The list of runs can be narrowed down to specific steps. All runs of the specified step will be submitted to the cluster. Also, individual runs IDs (step/run) can be used for submission.

  positional arguments:
    run                Submit only these runs to the cluster.

  optional arguments:
    -h, --help         show this help message and exit
    --even-if-dirty    This option must be set if the local git repository contains uncommited changes.
                       Otherwise uap will not run.
    --no-tool-checks   This option disables the otherwise mandatory checks for tool availability and version
    --cluster CLUSTER  Specify the cluster type. Default: [auto].
    --legacy           Use none array cluster submission.
    --force            Force to overwrite changed tasks.
    --ignore           Ignore chages of tasks and consider them finished.

.. _uap-fix-problems:

``fix-problems`` Subcommand
---------------------------

The ``fix-problems`` subcommand removes temporary files written by **uap** if
they are not required anymore.

Here is the usage information::

  $ uap <project-config>.yaml fix-problems -h
  usage: uap [<project-config>.yaml] fix-problems [-h] [--even-if-dirty]
                                                  [--no-tool-checks]
                                                  [--cluster CLUSTER]
                                                  [--first-error] [--details]
                                                  [--srsly]

  optional arguments:
    -h, --help         show this help message and exit
    --even-if-dirty    This option must be set if the local git repository contains uncommited changes.
                       Otherwise uap will not run.
    --no-tool-checks   This option disables the otherwise mandatory checks for tool availability and version
    --cluster CLUSTER  Specify the cluster type. Default: [auto].
    --first-error      Print stderr of the first failed cluster job.
    --details          Displays information about the files causing problems.
    --srsly            Delete problematic files.


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
  Hint: Run 'uap <project-config>.yaml fix-problems --first-error' to investigate what happended.
  Hint: Run 'uap <project-config>.yaml fix-problems --srsly' to fix these problems
        (that is, delete all problematic ping files).

Be nice and do as you've told.
Now you are able to resubmit your runs to the cluster.
You've fixed the problem, haven't you?

.. _uap-volatilize:

``volatilize`` Subcommand
-------------------------

The ``volatilize`` subcommand is useful to reduce the required disk space of
your analysis.
It works only if the :ref:`_volatile <config_file_volatile>` keyword is set in
the :ref:`analysis configuration file <analysis_configuration>` for.
As already mentioned there, steps marked as ``_volatile`` compute their output
files as normal but can be replaced by placeholder files if their dependent
steps are finished.

This subcommand provides usage information::

  $ uap <project-config>.yaml volatilize -h

  usage: uap [<project-config>.yaml] volatilize [-h] [--even-if-dirty]
                                                [--details] [--srsly]

  Save disk space by volatilizing intermediate results. Only steps marked with '_volatile: True' are considered.

  optional arguments:
    -h, --help       show this help message and exit
    --even-if-dirty  This option must be set if the local git repository
                     contains uncommited changes.
                     Otherwise uap will not run.
    --details        Shows which files can be volatilized.
    --srsly          Replaces files marked for volatilization with a placeholder.

After running ``volatilize --srsly`` the output files of the volatilized step
are replaced by placeholder files.
The placeholder files have the same name as the original files suffixed with
``.volatile.placeholder.yaml``.

.. _uap-render:

``render`` Subcommand
---------------------

The ``render`` subcommand generates graphs using graphviz.
The graphs either show the complete analysis or the execution of a single run.
At the moment ``--simple`` only has an effect in combination with ``--steps``.

This subcommand provides usage information::

  $ uap <project-config>.yaml render -h
  usage: uap [<project-config>.yaml] render [-h] [--even-if-dirty] [--files]
                                            [--steps] [--simple]
                                            [--orientation {left-to-right,right-to-left,top-to-bottom}]
                                            [run [run ...]]

  'render' generates DOT-graphs. Without arguments
  it takes the annotation file of each run and generates a graph,
  showing details of the computation.

  positional arguments:
    run                   Render only graphs for these runs.

  optional arguments:
    -h, --help            show this help message and exit
    --even-if-dirty       This option must be set if the local git repository
                          contains uncommited changes.
                          Otherwise uap will not run.
    --files               Renders a graph showing all files of the analysis.
                          [Not implemented yet!]
    --steps               Renders a graph showing all steps of the analysis and
                          their connections.
    --simple              Simplify rendered graphs.
    --orientation {left-to-right,right-to-left,top-to-bottom}
                          Defines orientation of the graph.
                          Default: 'top-to-bottom'

.. |argparse_link| raw:: html

   <a href="https://docs.python.org/2.7/library/argparse.html" target="_blank">argparse</a>

.. |sigkill_link| raw:: html

   <a href="https://en.wikipedia.org/wiki/Unix_signal#SIGKILL" target="_blank">SIGKILL</a>

.. |sigterm_link| raw:: html

   <a href="https://en.wikipedia.org/wiki/Unix_signal#SIGTERM" target="_blank">SIGTERM</a>
