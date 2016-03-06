..
  This is the documentation for uap. Please keep lines under 80 characters if
  you can and start each sentence on a new line as it decreases maintenance
  and makes diffs more readable.

.. title:: Command-Line Usage of uap

..
  This document aims to describe how to use **uap** via the command-line.

.. _cli_usage_uap

Command-Line Usage of **uap**
=============================

The **uap** software is meant to execute and monitor a data analysis which has
been configured in a YAML file (see :doc:`configuration`).
To start **uap** for a specific analysis.
The user needs to provide the path to a YAML configuration file as mandatory
parameter::

  $ uap <project-config>.yaml <subcommand> <options>


*<project-config>.yaml*
  Path to YAML file containing configuration of the data analysis

*<subcommand>*
  Any of the subcommands explained below (see `Explanation of Subcommands`_)

*<options>*
  Specific options for the used subcommand

Once the analysis is configured, several **uap** subcommands can be used
to interact with the pipeline. 
Everytime **uap** is started (with a subcommand) several things happen:

* The configuration file is read
* The tools given in the **Insert Link to Configuration:Tools Section**
  section are checked
* The input files are checked
* All tasks are calculated. 
 
If any of these steps fails, **uap** will print an error message with a
backtrace and it will crash.
This may seem a bit harsh, but after all, it's better to fail early than
to fail late if failing is unavoidable.

**uap** creates a symbolic link, if it does not exist yet, called 
``<configuration-file>-out``.
It points to the output directory defined in the <configuration-file>.
This symbolic link is located in the directory containing the
<configuration-file>.

There are a couple of global command line parameters which are valid for all 
scripts (well, actually, it's only one):

* ``--even-if-dirty``:
    Before doing anything else, the pipeline checks whether its source code 
    has been modified in any way via Git. 
    If yes, processing is stopped immediately unless this flag is specified.
    If you specify the flag, the fact that the repository was dirty will be 
    recorded in all annotations which are produces *including* a full Git diff.

..
    * ``--test-run``:
        When this parameter is specified, a ``head`` step is placed before all 
        first-level steps in the step tree, which returns the first 1000 lines 
        of every input file. 
        That way, a pipeline can be tested very quickly with a small input data 
        set.

The subcommands are described in detail below.

.. _ExplanationOfSubcommands:
Explanation of Subcommands
**************************

.. _uap_status

status
------

The status subcommand lists all tasks resulting from the configured steps and 
input samples.
At any time, each task is in one of the following states:

* **waiting** -- the task is waiting for input files to appear, or its input
  files are not up-to-date regarding their respective dependencies
* **ready** -- all input files are present and up-to-date regarding their 
  upstream input files (and so on, recursively), the task is ready and can 
  be started
* **queued** -- the task is currently queued and will be started "soon" 
  (if you use a computing cluster)
* **executing** -- the task is currently running on this or another machine
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



..
    Here is another example output with ``--test-run`` specified on the command 
    line. 
    Here, all top-level steps are prepended with a ``head`` step, which is 
    reflected in the task IDs::

        $ ./status.py --test-run
        [r] head/cutadapt/RIB0000784
        [r] head/cutadapt/RIB0000770
        [w] head/cutadapt/RIB0000784-R1
        [w] head/cutadapt/RIB0000784-R2
        [w] head/cutadapt/RIB0000770-R2
        [w] head/cutadapt/RIB0000770-R1
        [w] head/cutadapt/fix_cutadapt/RIB0000770
        [w] head/cutadapt/fix_cutadapt/RIB0000784
        tasks: 8 total, 2 ready, 6 waiting

Detailed information about a specific task can be obtained by specifying the 
task ID on the command line::

    $ uap <project-config>.yaml status cutadapt/Sample_COPD_2023-R1
    info:
      adapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG
    read_number: R1
    output_files:
      log:
        /home/michael/Desktop/rnaseq-pipeline/out/cutadapt-7708/Sample_COPD_2023-cutadapt-R1-log.txt:
        - /home/michael/Desktop/rnaseq-pipeline/copd-small/Sample_COPD_2023_R1.fastq.gz
      reads:
        /home/michael/Desktop/rnaseq-pipeline/out/cutadapt-7708/Sample_COPD_2023-cutadapt-R1.fastq.gz:
        - /home/michael/Desktop/rnaseq-pipeline/copd-small/Sample_COPD_2023_R1.fastq.gz
    state: FINISHED

This data structure is called the "run info" of a certain run and it 
represents a kind of plan which includes information about which output 
files will be generated and which input files they depend on -- this is 
stored in ``output_files``. 
Furthermore, necessary information for actually executing the task are 
recorded in ``info``. 
In this case, the final adapter has been determined by replacing ``((INDEX))`` 
in the configuration file's ``adapter-R1`` with the actual barcode index of 
the sample.

Because source steps produce no runs and therefore no tasks, they don't 
appear in the list produced by ``status.py``.
To see their task IDs, specify ``--sources``::

    $ uap <project-config>.yaml status --sources
    samples/Sample_COPD_2023
    
You can then specify the ID of a source task like the ID of any other task
to see its details::

    $ uap <project-config>.yaml status samples/Sample_COPD_2023
    info:
      index: ACAGTG
      paired_end: true
      read_number:
        Sample_COPD_2023_R1.fastq.gz: R1
        Sample_COPD_2023_R2.fastq.gz: R2
    output_files:
      reads:
        /home/michael/Desktop/rnaseq-pipeline/copd-small/Sample_COPD_2023_R1.fastq.gz: []
        /home/michael/Desktop/rnaseq-pipeline/copd-small/Sample_COPD_2023_R2.fastq.gz: []
      state: FINISHED



run-locally
-----------

The ``run-locally`` subcommand runs all non-finished tasks (or a subset) 
sequentially on the local machine. 
Feel free to cancel this script at any time, it won't put your project in a 
confused state.
However, if the ``run-locally`` subcommand receives a SIGKILL signal, the 
currently executing job will continue to run and the corresponding task
will be reported as ``executing`` by calling ``status`` subcommand for five more
minutes (SIGTERM should be fine and exit gracefully but *doesn't just yet*).
After that time, you will be warned that a job is marked as being currently
run but no activity has been seen for a while, along with further 
instructions about what to do in such a case (don't worry, it shouldn't 
happen by accident).

To execute one or more certain tasks, specify the task IDs on the command 
line. 
To execute all tasks of a certain step, specify the step name on the command 
line.

This subcommand provides usage information::
    
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
  The step_name is the name of an entry in the 'steps:' section as defined in '<project-config>.yaml'. A specific task is defined via its task ID 'step_name/run_id'. A list of all task IDs is returned by running:
  $ uap <project-config>.yaml status

  positional arguments:
    step_task        Can take multiple step names as input. A step name is the name of any entry in the 'steps:' section as defined in '<config>.yaml'. A list of all task IDs is returned by running:
                   $ uap <project-config>.yaml status.

  optional arguments:
    -h, --help       show this help message and exit
    --even-if-dirty  Must be set if the local git repository contains uncommited changtes. Otherwise the pipeline will not start.


.. NOTE:: Why is it safe to cancel the pipeline? 
    The pipeline is written in a way which expects processes to fail or 
    cluster jobs to disappear without notice. 
    This problem is mitigated by a design which relies on file presence and 
    file timestamps to determine whether a task is finished or not. 
    Output files are automatically written to temporary locations and later 
    moved to their real target directory, and it is not until the last file 
    rename operation has finished that a task is regarded as finished.
    
.. _uap_submit_to_cluster

submit-to-cluster
-----------------

The ``submit-to-cluster`` subcommand determines which tasks still have to be 
carried out and which supported cluster engine is available.
It then submits the jobs to the cluster if a cluster engine has been found. 
Dependencies are passed to cluster engine in a way that jobs that depend on
other jobs won't get scheduled until their dependencies have been satisfied. 
The files ``qsub-template.sh`` and ``sbatch-template.sh`` are used to submit
jobs, with ``#{ }`` fields being substituted with appropriate values.
Each submitted job calls **uap** with the ``run-locally`` subcommand on the
cluster nodes where the jobs are then run locally.

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
    
    $ ./run-locally.py -h
    usage: submit-to-cluster.py [-h] [--highmem] [--even-if-dirty]
                                [-s [STEP [STEP ...]]] [-t [TASK [TASK ...]]]

    This script submits all tasks configured in config.yaml to a Sun GridEngine 
    cluster via qsub. The list of tasks can be narrowed down by specifying a 
    step name (in which case all runs of this steps will be considered) or 
    individual tasks (step_name/run_id).

    optional arguments:
      -h, --help            show this help message and exit
      --highmem             this flag must be set if the highmem node of the 
                            cluster is being used.
      --even-if-dirty       Must be set if the local git repository contains 
                            uncommited changes. Otherwise the pipeline will not 
                            start.
      -s [STEP [STEP ...]], --step [STEP [STEP ...]]
                            Can take multiple step names as input. A step name 
                            is the name of any entry in the 'steps:' section as 
                            defined in 'config.yaml'
      -t [TASK [TASK ...]], --task [TASK [TASK ...]]
                            Can take multiple task ID(s) as input. A task ID 
                            looks like ths 'step_name/run_id'. A list of all 
                            task IDs is returned by running './status.py'.
fix-problems.py
---------------



render.py
---------


volatilize.py
-------------


