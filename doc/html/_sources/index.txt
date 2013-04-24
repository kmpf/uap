.. title:: rnaseq-pipeline

Introduction
============

The aim of this data processing pipeline is to enable simple and robust 
bioinformatics data evaluation.

**Simplicity:**

* The entire processing pipeline is described via a configuration file. Steps 
  are defined in a tree, and output files are written into a directory 
  structure mirroring this tree.
* Interaction with the pipeline happens through simple scripts which are used 
  to monitor the state of the pipeline and execute individual or all 
  remaining steps.
* To add a new processing step, a single Python file must be placed in 
  ``include/step`` which defines a class with two functions, one for 
  planning all jobs based on a list of input files or runs and possibly 
  additional information from previous steps and another function for 
  running a specific job.

**Robustness:**

* All steps write their output files to a temporary location. Only if a step 
  has completed successfully, the output files are copied to the correct 
  output directory.
* The output directory names are suffixed with a four-character hashtag 
  which mirrors the options specified for the step.
* Processing can be aborted and continued from the command line at any time. 
  This way, cluster failures are less critical because output files do not
  get compromised.
* Comprehensive annotations are written to the output directories, allowing 
  for later investigation about what exactly happened.
* Errors are caught as early as possible. Tools are checked for availability, 
  and the entire processing pipeline is calculated in advance before 
  jobs are being started or submitted to a cluster.

A pipeline is defined by two aspects:

* all processing steps arranges in a dependency tree
* its input samples

The combination of *steps* and *samples* result in a list of *tasks*, which 
can be executed sequentially or can be submitted to a cluster.

.. NOTE:: The design decision that steps are defined as a tree instead 
   of a full directed acyclic graph means that a step cannot have more than 
   one direct parent, like a directory in a file system cannot have more than 
   one parent directory. This means that a step cannot use the output of two 
   different steps as its input. However, any step may have more than one
   input or output file.

Setup
=====

The repository can be obtained like this::

    $ git clone spechtm@bioinf1:/home/spechtm/rnaseq-pipeline.git

After cloning the repository, run the bootstrapping script to create the 
required Python environment (which will be located in ``./python_env/``)::

    $ ./bootstrap.sh

There's no harm in accidentally running this script multiple times. Also,
it will compile ``cat4m``, a tool which can be found at ``./tools/cat4m``
and which is able to read arbitrary input files in chunks of 4 MB
and prints them to stdout.

The configuration file
----------------------

Next, edit ``config.sample.yaml`` and save it as ``config.yaml``. Although 
writing the configuration may seem a bit complicated, the trouble pays off 
later because further interaction with the pipeline is quite simple. Here is 
a sample configuration:

.. code-block:: yaml

    # This is the rnaseq-pipeline configuration file.
    email: micha.specht@gmail.com
    sources:
    - run_folder_source: { path: in }
    destination_path: out
    steps: |
        - cutadapt {
            adapter-R1: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC((INDEX))ATCTCGTATGCCGTCTTCTGCTTG"
            adapter-R2: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
        }
            - fix_cutadapt
    tools:
        cutadapt:
            path: 'tools/cutadapt-1.2.1/bin/cutadapt'
            get_version: '--version'
        pigz:
            path: 'pigz'
            get_version: '--version'
        dd:
            path: 'dd'
            get_version: '--version'

In the configuration, the following aspects of the pipeline are defined:

* ``sources`` -- there can be multiple sources of different types:

  * run folders
  * plain fastq.gz files with additional information
  
* ``destination_path`` -- this is where result files, annotations and 
  temporary files are written to
* ``steps`` -- defines the processing step arranged in a tree
* ``tools`` -- defines all tools used in the pipeline and how to determine 
  their versions (for later reference)
* ``email`` -- when submitting jobs on a cluster, messages will be sent to 
  this email address (nobody@example.com by default)
  
Scripts
=======

Once the project is set up, there are several scripts which can be used to 
execute and monitor the pipeline. All scripts have a couple of properties in 
common:

* On startup, the configuration is read, tools are checked, input files are 
  collected, and all tasks are calculated. If any of these steps fails, the 
  script will print an error message with a backtrace and it will crash.
* For convenience, a symbolic link called ``out`` will be placed in the 
  pipeline's directory which points to the output directory defined in the 
  configuration file. If ``out`` already exists, it is left untouched.

There are a couple of global command line parameters which are valid for all 
scripts:

* ``--even-if-dirty``:
    Before doing anything else, the pipeline checks whether its source code 
    has been modified in any way via Git. If yes, processing is stopped 
    immediately unless this flag is specified.
* ``--test-run``:
    When this parameter is specified, a ``head`` step is placed before all 
    first-level steps in the step tree, which returns the first 1000 lines 
    of every input file. That way, a pipeline can be tested very quickly 
    with a small input data set.

In the following, the scripts are described in detail.

status.py
---------

The status script lists all tasks resulting from the configured steps and 
input samples. At the beginning of each line, the status of each task is 
denoted by [w], [r], and [f], corresponding to:

* **waiting** -- the taks is waiting for input files to appear or to be updated
* **ready** -- all input files are present and up-to-date regarding their 
  upstream input files, task can be started
* **finished** -- all output files are in place and up-to-date

.. NOTE:: In the current design, there is no mechanism to indicate 
    whether a task is currently running or has been submitted to a cluster.

Here is an example output::

    $ ./status.py
    [r] cutadapt/RIB0000784-R1
    [r] cutadapt/RIB0000784-R2
    [r] cutadapt/RIB0000770-R2
    [r] cutadapt/RIB0000770-R1
    [w] cutadapt/fix_cutadapt/RIB0000770
    [w] cutadapt/fix_cutadapt/RIB0000784
    tasks: 6 total, 4 ready, 2 waiting

Here is another example output with ``--test-run`` specified on the command 
line. Here, all top-level steps are prepended with a ``head`` step, which is 
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

    $ ./status.py cutadapt/RIB0000770-R1
    info:
    adapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG
    read: R1
    output_files:
    log:
        out/cutadapt-7708988d/RIB0000770-cutadapt-R1-log.txt:
        - in/Unaligned/Project_A/Sample_RIB0000770/RIB0000770_TAGCTT_L001_R1_001.fastq.gz
        - in/Unaligned/Project_A/Sample_RIB0000770/RIB0000770_TAGCTT_L002_R1_001.fastq.gz
        - in/Unaligned/Project_A/Sample_RIB0000770/RIB0000770_TAGCTT_L003_R1_001.fastq.gz
        - in/Unaligned/Project_A/Sample_RIB0000770/RIB0000770_TAGCTT_L004_R1_001.fastq.gz
        - in/Unaligned/Project_A/Sample_RIB0000770/RIB0000770_TAGCTT_L005_R1_001.fastq.gz
        - in/Unaligned/Project_A/Sample_RIB0000770/RIB0000770_TAGCTT_L006_R1_001.fastq.gz
        - in/Unaligned/Project_A/Sample_RIB0000770/RIB0000770_TAGCTT_L007_R1_001.fastq.gz
        - in/Unaligned/Project_A/Sample_RIB0000770/RIB0000770_TAGCTT_L008_R1_001.fastq.gz
    reads:
        out/cutadapt-7708988d/RIB0000770-cutadapt-R1.fastq.gz:
        - in/Unaligned/Project_A/Sample_RIB0000770/RIB0000770_TAGCTT_L001_R1_001.fastq.gz
        - in/Unaligned/Project_A/Sample_RIB0000770/RIB0000770_TAGCTT_L002_R1_001.fastq.gz
        - in/Unaligned/Project_A/Sample_RIB0000770/RIB0000770_TAGCTT_L003_R1_001.fastq.gz
        - in/Unaligned/Project_A/Sample_RIB0000770/RIB0000770_TAGCTT_L004_R1_001.fastq.gz
        - in/Unaligned/Project_A/Sample_RIB0000770/RIB0000770_TAGCTT_L005_R1_001.fastq.gz
        - in/Unaligned/Project_A/Sample_RIB0000770/RIB0000770_TAGCTT_L006_R1_001.fastq.gz
        - in/Unaligned/Project_A/Sample_RIB0000770/RIB0000770_TAGCTT_L007_R1_001.fastq.gz
        - in/Unaligned/Project_A/Sample_RIB0000770/RIB0000770_TAGCTT_L008_R1_001.fastq.gz

The details of this data structure are explained below. It represents a kind 
of plan which includes information about which output files will be generated 
and which input files they depend on -- this is stored in ``output_files``. 
Furthermore, necessary information for actually executing the task are 
recorded in ``info``. In this case, the final adapter has been determined by 
replacing ``((INDEX))`` in the configuration file's ``adapter-R1`` with the 
actual barcode index of the sample.


run-locally.py
--------------

The ``run-locally.py`` script runs all non-finished tasks (or a subset) 
sequentially on the local machine. Feel free to cancel this script at any 
time, it won't put your project in a confused state.

To execute one or more certain tasks, specify the task IDs on the command 
line.

.. NOTE:: Why is it safe to cancel the pipeline? The pipeline is 
    written in a way which expects processes to fail or cluster jobs to 
    disappear without notice. This problem is mitigated by a design which 
    relies on file presence and file timestamps to determine whether a task 
    is finished or not. Output files are automatically written to temporary 
    locations and later moved to their real target directory, and it is not 
    until the last file rename operation has finished that a task is 
    regarded as finished.

submit-to-cluster.py
--------------------

The ``submit-to-cluster.py`` script determines which tasks still have to be 
carried out and submits the jobs to a GridEngine cluster by calling ``qsub``. 
Dependencies are passed to ``qsub`` via the ``-hold_jid`` option, which means 
that jobs that depend on other jobs won't get scheduled until their 
dependencies have been satisfied. The file ``qsub-template.sh`` is used to 
submit jobs, with ``#{ }`` fields being substituted with appropriate values.

The file ``quotas.yaml`` can be used to define different quotas for different 
systems:

.. code-block:: yaml

    "frontend[12]":
        default: 5

In the example above, a default quota of 5 is defined for hosts with a 
hostname of ``frontend1`` or ``frontend2`` (the name is a regular expression). 
Different quotas can be defined for each step. A quota of 5 means that no 
more than 5 jobs of on kind (the same step) will be run in parallel.

  
Sources
=======

In the following, the different types of sources are described in detail.

Run folder source
-----------------

Here's an example:

.. code-block:: yaml

    - run_folder_source: { path: in }

This source looks for fastq.gz files in ``[path]/Unaligned/Project_*/Sample_*`` 
and pulls additional information from CSV sample sheets it finds. It also 
makes sure that index information for all samples is coherent and unambiguous.

FASTQ source
------------

Here's an example:

.. code-block:: yaml

    - fastq_source:
        pattern: /data/original-fastq/&#42;.fastq.gz
        group: (Sample_COPD_\d+)_R[12].fastq.gz
        indices: copd-barcodes.csv

Input files are collected as defined by ``pattern`` and grouped into samples 
according to ``group``, which is a regular expression. All groups defined in 
the regex ``(  )`` are used to construct the sample name, here it is used to 
declare that both R1 and R2 files belong to the same sample. Indices are 
read from the CSV file specified by ``indices``.

..    
    .. automodule:: abstract_source

    .. autoclass:: AbstractSource
        :members:
    
Steps
=====

Steps are defined in a dependency tree. However, the syntax is a bit peculiar: 
The ``|`` after ``steps:`` is YAML-specific syntax and it defines a string 
spanning multiple lines in which line breaks and indentation is maintained. 
The string is later parsed by the pipeline and the most important parts are 
the individual steps which are to be performed. The relationship betweens 
steps is declared via indentation.

.. NOTE:: Why do we need the ``|`` symbol in the steps definition? Neither the 
    list nor the dictionary syntax allow for a concise definition of a step 
    tree with options. Think of the step definition as a nested list with an 
    option hash attached to every item.

Steps may have options, which must be placed in between ``{`` curly braces 
``}``. Options can be specified on a single line (in this case, individual 
key/value pairs must be separated by comma) or may span multiple lines, 
following standard YAML block syntax.

..
    .. automodule:: abstract_step
        :members:

    .. autoclass:: AbstractStep
        :members:
    
Miscellaneous
-------------

Head
~~~~
    
.. graphviz::

    digraph foo {
        rankdir=LR;
        splines=true;
        graph [fontname = Helvetica, fontsize = 12, nodesep = 0.2, ranksep = 0.3];
        node [fontname = Helvetica, fontsize = 12, shape = rect, style=filled, color="#404040", fillcolor="#ffffff"];
        edge [fontname = Helvetica, fontsize = 12, color="#404040"];

        head [fillcolor = "#fce94f", color = "#c4a000"];
        in_any [label = "any (*)"];
        out_any [label = "any (*)"];
        
        in_any -> head;
        head -> out_any;
        
    }

.. autosimpleclass:: head.Head
    
Break
~~~~~

.. graphviz::

    digraph foo {
        rankdir=LR;
        splines=true;
        graph [fontname = Helvetica, fontsize = 12, nodesep = 0.2, ranksep = 0.3];
        node [fontname = Helvetica, fontsize = 12, shape = rect, style=filled, color="#404040", fillcolor="#ffffff"];
        edge [fontname = Helvetica, fontsize = 12, color="#404040"];

        break [fillcolor = "#fce94f", color = "#c4a000"];
        in_any [label = "any (*)"];
        
        in_any -> break;
    }

.. autosimpleclass:: break.Break
    
Preprocessing
-------------
    
Adapter clipping
~~~~~~~~~~~~~~~~

Cutadapt
^^^^^^^^

.. graphviz::

    digraph foo {
        rankdir=LR;
        splines=true;
        graph [fontname = Helvetica, fontsize = 12, nodesep = 0.2, ranksep = 0.3];
        node [fontname = Helvetica, fontsize = 12, shape = rect, style=filled, color="#404040", fillcolor="#ffffff"];
        edge [fontname = Helvetica, fontsize = 12, color="#404040"];

        cutadapt [fillcolor = "#fce94f", color = "#c4a000"];
        in_reads [label = "reads\n(fastq.gz)"];
        out_reads [label = "reads\n(fastq.gz)"];
        out_log [label = "log\n(.txt)"];
        
        in_reads -> cutadapt;
        cutadapt -> out_reads;
        cutadapt -> out_log;
    }

.. autosimpleclass:: cutadapt.Cutadapt
    
Fix cutadapt
^^^^^^^^^^^^

.. graphviz::

    digraph foo {
        rankdir=LR;
        splines=true;
        graph [fontname = Helvetica, fontsize = 12, nodesep = 0.2, ranksep = 0.3];
        node [fontname = Helvetica, fontsize = 12, shape = rect, style=filled, color="#404040", fillcolor="#ffffff"];
        edge [fontname = Helvetica, fontsize = 12, color="#404040"];

        fix_cutadapt [fillcolor = "#fce94f", color = "#c4a000"];
        in_reads [label = "reads\n(fastq.gz)"];
        out_reads [label = "reads\n(fastq.gz)"];
        
        in_reads -> fix_cutadapt;
        fix_cutadapt -> out_reads;
    }

.. autosimpleclass:: fix_cutadapt.FixCutadapt
    
Aligners
--------

Segemehl
~~~~~~~~

.. graphviz::

    digraph foo {
        rankdir=LR;
        splines=true;
        graph [fontname = Helvetica, fontsize = 12, nodesep = 0.2, ranksep = 0.3];
        node [fontname = Helvetica, fontsize = 12, shape = rect, style=filled, color="#404040", fillcolor="#ffffff"];
        edge [fontname = Helvetica, fontsize = 12, color="#404040"];

        segemehl [fillcolor = "#fce94f", color = "#c4a000"];
        in_reads [label = "reads\n(fastq.gz)"];
        out_mapped_reads [label = "mapped reads\n(sam.gz)"];
        out_unmapped_reads [label = "unmapped reads\n(fastq.gz)"];
        out_log [label = "log\n(txt)"];
        
        in_reads -> segemehl;
        segemehl -> out_mapped_reads;
        segemehl -> out_unmapped_reads;
        segemehl -> out_log;
    }

.. autosimpleclass:: segemehl.Segemehl

Extending the pipeline
======================

Checklist
---------

Here's a couple of things which should be kept in mind when implementing new 
steps or modifying existing steps:

* Make sure errors already show up in ``setup_runs`` instead of ``execute``.
  That way, wasting precious cluster waiting time is avoided. Look out for 
  things that may fail, and do them in ``setup_runs``. Use the ``info`` 
  entry in the returned ``run_info`` structure to pass the resulting 
  information to ``execute``.
* Likewise, make sure that the tools you'll need in execute are already 
  available in ``setup_runs``::
  
    # make sure tools are available
    self.tool('pigz')
    self.tool('cutadapt')
    
* Make sure your disk access is as cluster-friendly as possible (which 
  primarily means using large block sizes and preferably no seek operations). 
  If possible, use ``unix_pipeline`` to wrap your commands in ``pigz``, 
  ``dd``, or ``cat4m`` with a block size of at least 4 MB. Although this is 
  not possible in every case (for example when seeking inside files is 
  involved), it is straightforward with tools that read a continuous stream 
  from ``stdin`` and write a continuous stream to ``stdout``.

To-do list
==========

Timestamps:
    ``unix_pipeline`` log messages should include timestamps.

Getting started package:
    We need a small package which demonstrates a quick pipeline, including
    the configuration and all required tools.
    
Capture process output:
    For all processes launched via ``unix_pipeline``, the respective stdout 
    and stderr should be recorded and the last kB remembered, so that it can
    be included into the error message if a pipeline fails. Also, the 
    captured output should be incorporated into the YAML annotations which
    are written for every output file.
    
    *Plus:* This would also allow for the automatic generation of SHA1 
    checksums on the fly.
    
Steps should be able to access all ancestors:
    All upstream steps should be accessible via their step name or output 
    file key.
    
On-the-fly steps:
    We need a way to skip writing certain output files and have them flow 
    temporarily through a pipe only, if possible. This is a disk space-saving 
    feature only and has no effect on the outcome of the pipeline. However,
    it would require that a step is capable of being run *on-the-fly* which
    means it must read and write in a single stream.
    
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

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. toctree::
    :maxdepth: 2

    api