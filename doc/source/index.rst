..
  This is the documentation for rnaseq-pipeline. Please keep lines under
  80 characters if you can and start each sentence on a new line as it 
  decreases maintenance and makes diffs more readable.
  
.. title:: rnaseq-pipeline


Introduction
============

The **rnaseq-pipeline** package is a framework to configure, run, and control
high-throughput sequencing analyses.

The aim of this data processing pipeline is to enable robust and straightforward
bioinformatics data evaluation. It is implemented in Python, runs under
GNU/Linux and can be controlled from the command-line interface. Although the
primary focus is the evaluation of RNASeq data, its design allows for a variety
of other applications.

**Table of contents**

.. toctree::
   :maxdepth: 2

   documentation
   steps
   api


General usage
-------------

This package *does not* provide a number of tools which are downloaded and
installed system-wide to provide certain functioniality.
The intention of this system is to provide a robust and traceable framework
for data evaluation in scientific experiments which uses other tools and
manages individual data processing steps and their inter-dependencies.
    
The recommended workflow for running a data evaluation for an experiment is 
as follows:

1. Check-out the rnaseq-pipeline repository via Git.
2. Setup the project by writing the configuration file.
3. Add steps or other functionality as needed (optional).
4. Run the pipeline.
5. Have your changes (if there are any) merged back into the main repository,
   to the advantage of the scientific community.

This leaves you with:

* Your original input files, which are left untouched.
* The experiment-specific pipeline repository.  
  You should keep this repository for later reference and you could even
  make it publicly available along with your input files for anybody to
  re-run the entire data evaluation or parts thereof.
* The output directory containing all output files and comprehensive 
  annotations.
  These annotations include detailed information for every output file,
  including which steps have been executed and the Git SHA1 hash of
  the pipeline repository at the time the data processing took place.
  In many cases, these annotations also include information about all
  inter-process streams and output files, including SHA1 checksums, file 
  sizes, and line counts.

Core aspects
------------

**Robustness:**

* All steps write their output files to a temporary location. 
  Only if a step has completed successfully, the output files are copied to 
  the correct output directory.
* The output directory names are suffixed with a four-character hashtag 
  which mirrors the options specified for the step.
* Processing can be aborted and continued from the command line at any time. 
  This way, cluster failures are less critical because output files do not
  get compromised.
* Errors are caught as early as possible. Tools are checked for availability, 
  and the entire processing pipeline is calculated in advance before 
  jobs are being started or submitted to a cluster.
  
**Traceability:**

* Comprehensive annotations are written to the output directories, allowing 
  for later investigation about what exactly happened.
      
**Simplicity:**

* The entire processing pipeline is described via a configuration file. 
  Steps are defined in a directed acyclic graph (DAG).
* Interaction with the pipeline happens through a handful of scripts which 
  are used to monitor the state of the pipeline and execute individual or all 
  remaining steps.

Design
------

The central part of the pipeline is its definition of the steps which are to 
be carried out.
Steps are organized in a dependency graph (a directed acyclic graph) -- every 
step may have one or more parent steps, which may in turn have other parent 
steps, and so on.
Steps without parents are usually sources which provide source files, for
example FASTQ files with the raw sequences obtained from the sequencer,
genome sequence databases or annotation tracks.

Each step defines a number of runs and each run represents a piece of the
entire data evaluation, typically at the level of a single sample.
A certain *run* of a certain *step* is called a *task*.
While the steps only describe what needs to be done on a very abstract level,
it is through the individual runs of each step that a pipeline-wide list of 
actual tasks becomes available.
Each run may provide a number of output files which depend on output files
of one or several runs from parent steps.

Source steps define a run for every input sample, and a subsequent step
may:

* define the same number of runs, 
* define more runs (for example when R1 and R2 reads in a paired-end RNASeq 
  experiment should be treated separately),
* define fewer runs (usually towards the end of a pipeline, where results are
  summarized).

=======
Remarks
=======

This documentation has been created using `Sphinx <http://sphinx-doc.org/>`_
and `reStructuredText <http://docutils.sourceforge.net/rst.html>`_.
=======
Setup
=====

The repository can be obtained like this::

    $ git clone git@github.com:tiennes/rnaseq-pipeline.git

After cloning the repository, run the bootstrapping script to create the 
required Python environment (which will be located in ``./python_env/``)::

    $ ./bootstrap.sh

There's no harm in accidentally running this script multiple times. 
Also, it will compile ``cat4m``, a tool which can be found at 
``./tools/cat4m`` and which is able to read arbitrary input files in chunks 
of 4 MB and print them to stdout (we'll need this often in the pipeline,
as ``cat`` reads in system-default blocks of 32 kB which is ok for a normal
system but leads to high I/O load on a cluster system).

The configuration file
----------------------

Next, edit ``config.sample.yaml`` and save it as ``config.yaml``. 
Although writing the configuration may seem a bit complicated, the trouble 
pays off later because further interaction with the pipeline is quite simple. 
Here is a sample configuration:

.. code-block:: yaml

    # This is the rnaseq-pipeline configuration file.

    destination_path: "/home/michael/test-pipeline/out"

    steps:
        fastq_source:
            pattern: /home/michael/test-pipeline/fastq/*.fastq.gz
            group: (Sample_COPD_\d+)_R[12]-head.fastq.gz
            indices: copd-barcodes.csv
            paired_end: yes
            
        cutadapt:
            _depends: fastq_source
            adapter-R1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC((INDEX))ATCTCGTATGCCGTCTTCTGCTTG
            adapter-R2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
            
    tools:
        cutadapt:
            path: /home/michael/Desktop/rnaseq-pipeline/tools/cutadapt-1.2.1/bin/cutadapt
            get_version: '--version'
        pigz:
            path: pigz
            get_version: '--version'
        head:
            path: head
            get_version: '--version'
        cat4m:
            path: ./tools/cat4m

In the configuration, the following aspects of the pipeline are defined:

* ``destination_path`` -- this is where result files, annotations and 
  temporary files are written to
* ``steps`` -- defines the processing step arranged in a DAG
* ``tools`` -- defines all tools used in the pipeline and how to determine 
  their versions (for later reference)
* ``email`` -- when submitting jobs on a cluster, messages will be sent to 
  this email address by the cluster scheduler (nobody@example.com by default)
  
Steps
~~~~~
  
Steps are defined in a directed acyclic graph. 
In the configuration, the ``steps`` dictionary contains a key for every
step, therefore each step must have a unique name.
There are two ways to name a step:

.. code-block:: yaml

    steps:
        # here, the step name is unchanged, it's a cutadapt step which is also called 'cutadapt'
        cutadapt:
            ... # options following
            
        # here, we also insert a cutadapt step, but we give it a different name: 'clip_adapters'
        clip_adapters (cutadapt):
            ... # options following
            
Source steps are special in the way that they provide files without doing
anything, and they are usually the first steps in a pipeline because they
have no dependencies.
Regular steps, on the other hand, need to define their dependencies via
the ``_depends`` key which may either be ``null``, a step name, or a list
of step names.

.. code-block:: yaml

    steps:
        # the source step which depends on nothing
        fastq_source:
            # ...
            
        # the first processing step, which depends on the sources
        cutadapt:
            _depends: fastq_source
        
        # the second processing step, which depends on the cutadapt step
        fix_cutadapt:
            _depends: cutadapt
                
If you want to cut off entire branches of the step graph, set the ``_BREAK`` 
flag in a step definition, which will force the step to produce no runs
(which will in turn give all following steps nothing to do, thereby 
effectively disabling these steps):
        

.. code-block:: yaml

    steps:
        fastq_source:
            # ...
            
        cutadapt:
            _depends: fastq_source
        
        # this step and all following steps will not be executed
        fix_cutadapt:
            _depends: cutadapt
            _BREAK: true

Tools
~~~~~

All tools which are used in the pipeline must be specified in the configuration 
file. The pipeline determines and records their versions for future reference.

By default, version determination is simply attempted by calling the program
without command-line arguments.

If a certain argument is required, specify it in ``get_version``. 
If a tool does not exit with exit code 0, find out which code it is by typing
``echo $?`` into Bash and specify the exit code in ``exit_code``.

.. code-block:: yaml

    tools:
        cutadapt:
            path: /home/michael/Desktop/rnaseq-pipeline/tools/cutadapt-1.2.1/bin/cutadapt
            get_version: '--version'
            
        head:
            path: head
            get_version: '--version'

Scripts
=======

Once the project is set up, there are several scripts which can be used to 
execute and monitor the pipeline. 
All scripts have a couple of properties in common:

* On startup, the configuration is read, tools are checked, input files are 
  collected, and all tasks are calculated. 
  If any of these steps fail, the script will print an error message with 
  a backtrace and it will crash.
  This may seem a bit harsh, but after all, it's better to fail early than
  to fail late if failing is unavoidable.
* For convenience, a symbolic link called ``out`` will be placed in the 
  pipeline's directory which points to the output directory defined in the 
  configuration file. 
  If ``out`` already exists, it is left untouched.

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

In the following, the scripts are described in detail.

status.py
---------

The status script lists all tasks resulting from the configured steps and 
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

    $ ./status.py
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

    $ ./status.py --summarize
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

    $ ./status.py --graph
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

    $ ./status.py cutadapt/Sample_COPD_2023-R1
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

    $ ./status.py --sources
    samples/Sample_COPD_2023
    
You can then specify the ID of a source task like the ID of any other task
to see its details::

    $ ./status.py samples/Sample_COPD_2023
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



run-locally.py
--------------

The ``run-locally.py`` script runs all non-finished tasks (or a subset) 
sequentially on the local machine. 
Feel free to cancel this script at any time, it won't put your project in a 
confused state.
However, if the ``run-locally.py`` script receives a SIGKILL signal, the 
currently executing job will continue to run and the corresponding task
will be reported as ``executing`` by ``status.py`` for five more minutes
(SIGTERM should be fine and exit gracefully but *doesn't just yet*).
After that time, you will be warned that a job is marked as being currently
run but no activity has been seen for a while, along with further 
instructions about what to do in such a case (don't worry, it shouldn't 
happen by accident).

To execute one or more certain tasks, specify the task IDs on the command 
line. 
To execute all tasks of a certain step, specify the step name on the command 
line.

This script provides usage information::
    
    $ ./run-locally.py -h

    usage: run-locally.py [-h] [--even-if-dirty] [-s [STEP [STEP ...]]]
                          [-t [TASK [TASK ...]]]

    This script starts the 'rnaseq-pipeline' on the local machine. It can be 
    used to start:
     * all tasks of the pipeline as configured in 'config.yaml'
     * all tasks defined by a specific step in 'config.yaml'
     * one or more steps

    To start the complete pipeline as configured in 'config.yaml' execute:
    $ ./run-locally.py

    To start a specific step execute:
    $ ./run-locally.py <step_name>

    To start a specific task execute:
    $ ./run-locally.py <step_name/run_id>

    The step_name is the name of an entry in the 'steps:' section as defined in 
    'config.yaml'. A specific task is defined via its task ID 'step_name/run_id'.
    A list of all task IDs is returned by running './status.py'.

    optional arguments:
      -h, --help            show this help message and exit
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


.. NOTE:: Why is it safe to cancel the pipeline? 
    The pipeline is written in a way which expects processes to fail or 
    cluster jobs to disappear without notice. 
    This problem is mitigated by a design which relies on file presence and 
    file timestamps to determine whether a task is finished or not. 
    Output files are automatically written to temporary locations and later 
    moved to their real target directory, and it is not until the last file 
    rename operation has finished that a task is regarded as finished.
    
submit-to-cluster.py
--------------------

The ``submit-to-cluster.py`` script determines which tasks still have to be 
carried out and submits the jobs to a GridEngine cluster via ``qsub``. 
Dependencies are passed to ``qsub`` via the ``-hold_jid`` option, which means 
that jobs that depend on other jobs won't get scheduled until their 
dependencies have been satisfied. 
The file ``qsub-template.sh`` is used to submit jobs, ``#{ }`` fields 
being substituted with appropriate values. Each submitted job calls the 
``run-locally.py`` script that actually tells the pipeline what to do.

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


This script provides usage information::
    
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





Annotations
===========
    
Upon successful completion of a task, an extensive YAML-formatted annotation 
is placed next to the output files in a file called 
``.[task_id]-annotation.yaml``.
Also, for every output file, a symbolic link to this file is created:
``.[output_filename].annotation.yaml``.

Finally, the annotation is rendered via GraphViz, if available.
Rendering can also be done at a later time using annotations as input.
The annotation can be used to determine at a later time what exactly happened.
Also, annotations may help to identify bottlenecks.

+---------------------------------------+-----------------------------------------------+
| .. image:: _static/cutadapt.png       | .. image:: _static/cpu-starving.png           |
|   :height: 500                        |   :height: 500                                |
|                                       |                                               |
| Annotation graph of a ``cutadapt``    | In this graph, it becomes evident that        |
| run. CPU and RAM usage for individual | the ``fix_cutadapt.py`` process in the middle |
| processes are shown, file sizes       | gets throttled by the following two ``pigz``  |
| and line counts are shown for         | processes, which only run with one core       |
| output files and inter-process        | each and therefore cannot compress the        |
| streams.                              | results fast enough.                          |
+---------------------------------------+-----------------------------------------------+
    
Steps
=====

A detailed description of all availble steps follows.
               
Source steps
------------

Source steps provide input files for the pipeline, such as RNA sequences.

Run folder source
~~~~~~~~~~~~~~~~~

Here's an example:

.. code-block:: yaml

    - run_folder_source: { path: in }

This source looks for fastq.gz files in 
``[path]/Unaligned/Project_*/Sample_*`` and pulls additional information from 
CSV sample sheets it finds. 
It also makes sure that index information for all samples is coherent and 
unambiguous.

FASTQ source
~~~~~~~~~~~~

Here's an example:

.. code-block:: yaml

    - fastq_source:
        pattern: /data/original-fastq/&#42;.fastq.gz
        group: (Sample_COPD_\d+)_R[12].fastq.gz
        indices: copd-barcodes.csv

Input files are collected as defined by ``pattern`` and grouped into samples 
according to ``group``, which is a regular expression. 
All groups defined in the regex ``(  )`` are used to construct the sample 
name, here it is used to declare that both R1 and R2 files belong to the 
same sample. 
Indices are read from the CSV file specified by ``indices``.

..    
    .. automodule:: abstract_source

    .. autoclass:: AbstractSource
        :members:

*(detailed step descriptions to follow...)*

..
    Miscellaneous
    -------------

    Head
    ~~~~
        
    .. autosimpleclass:: head.Head
        
    Preprocessing
    -------------
        
    Adapter clipping
    ~~~~~~~~~~~~~~~~

    Cutadapt
    ^^^^^^^^

    .. autosimpleclass:: cutadapt.Cutadapt
        
    Fix cutadapt
    ^^^^^^^^^^^^

    .. autosimpleclass:: fix_cutadapt.FixCutadapt
        
    Aligners
    --------

    Segemehl
    ~~~~~~~~

    .. autosimpleclass:: segemehl.Segemehl

Extending the pipeline
======================

To add a new processing step, a single Python file must be placed in 
``include/step`` which defines a class with two functions, one for planning 
all jobs based on a list of input files or runs and possibly additional 
information from previous steps and another function for running a specific 
job.


Checklist
---------

Here's a couple of things which should be kept in mind when implementing new 
steps or modifying existing steps:

* Make sure errors already show up in ``setup_runs`` instead of ``execute``.
  That way, wasting precious cluster waiting time is avoided. 
  Look out for things that may fail, and do them in ``setup_runs``. 
  Use the ``info`` entry in the returned ``run_info`` structure to pass the 
  resulting information to ``execute``.
* Likewise, make sure that the tools you'll need in execute are already 
  available in ``setup_runs``::
  
    # make sure tools are available
    self.tool('pigz')
    self.tool('cutadapt')
    
* Make sure your disk access is as cluster-friendly as possible (which 
  primarily means using large block sizes and preferably no seek operations). 
  If possible, use ``unix_pipeline`` to wrap your commands in ``pigz``, 
  ``dd``, or ``cat4m`` with a large block size like 4 MB. 
  Although this is not possible in every case (for example when seeking 
  in files is involved), it is straightforward with tools that read a 
  continuous stream from ``stdin`` and write a continuous stream to 
  ``stdout``.

To-do list
==========

Timestamps:
    ``unix_pipeline`` log messages should include timestamps.

Getting started package:
    We need a small package which demonstrates a quick pipeline, including
    the configuration and all required tools.
    
Steps should be able to access all ancestors:
    All upstream steps should be accessible via their step name or output 
    file key.
    
On-the-fly steps:
    We need a way to skip writing certain output files and have them flow 
    temporarily through a pipe only, if possible. 
    This is a disk space-saving feature only and has no effect on the 
    outcome of the pipeline. However, it would require that a step is 
    capable of being run *on-the-fly* which means it must read and write in 
    a single stream.
    
    Here's an example:
    
    .. graphviz::
        digraph foo {
            rankdir=LR;
            splines=true;
            graph [fontname = Helvetica, fontsize = 12, nodesep = 0.2, ranksep = 0.3];
            node [fontname = Helvetica, fontsize = 12, shape = rect, style=filled, color="#404040", fillcolor="#ffffff"];
            edge [fontname = Helvetica, fontsize = 12, color="#404040"];

            segemehl [fillcolor = "#fce94f", color = "#c4a000"];
            in_reads [label = "reads\n(fastq.gz)"];
            mapped_reads [label = "mapped reads\n(sam.gz)"];
            some_filter [fillcolor = "#fce94f", color = "#c4a000"];
            filtered_reads [label = "filtered reads\n(sam.gz)"];
            htseq_count [label = "htseq-count", fillcolor = "#fce94f", color = "#c4a000"];
            counts [label = "counts"];
            
            in_reads -> segemehl -> mapped_reads -> some_filter -> filtered_reads;
            filtered_reads -> htseq_count -> counts;

            subgraph cluster_food {
                some_filter; filtered_reads;
                label = "on-the-fly step, filtered reads\nnever get written to disk";
                graph [style=dashed, color="#808080"];
            }
        }
        
Miscellaneous input files:
    Genome files and their index such as used by segemehl should not be defined
    via a fixed path.
    For traceability, it would be preferable to specify the hg19.fa URL and
    checksum and have the index generated by a step which the segemehl step
    depends on.
    
Make ``run-locally.py`` exit gracefully on receiving SIGTERM.

Show statistics for executing tasks:
    When showing currently executing tasks, show how long this job has already been
    running and how it relates to jobs that have already finished.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
>>>>>>> working
