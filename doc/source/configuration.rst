..
  This is the documentation for rnaseq-pipeline. Please keep lines under
  80 characters if you can and start each sentence on a new line as it 
  decreases maintenance and makes diffs more readable.

.. title:: Documentation

..
  This document aims to describe how to configure **uap**.

Configuration of **uap**
========================

**uap** is made to control the execution of data analyses which are defined
in `YAML <http://www.yaml.org/>`_ files.
Each file describes a complete analysis.
Further on these files are called analysis file(s).

The analysis files consist of four sections (let's just call them sections,
although technically, they are keys):

* ``destination_path`` -- points to the directory where the result files,
  annotations and temporary files are written to
* ``email`` -- when submitting jobs on a cluster, messages will be sent to 
  this email address by the cluster scheduler (nobody@example.com by default)
* ``constants`` -- defines constants for later use (define repeatedly used
  values as constants to increase readability of the following sections)
* ``tools`` -- defines all tools used in the pipeline and how to determine 
  their versions (for later reference)
* ``steps`` -- defines the processing step and their order 

If you want to know more about the notation that is used in this file, have a
closer look at the `YAML definition <http://www.yaml.org/>`_.

Sections of config.yaml
***********************

destination_path
~~~~~~~~~~~~~~~~

The value of ``destination_path`` is the directory where **uap** is going
to store the created files. It is possible to use a different directory for
volatile files (see ).

.. code-block:: yaml

    destination_path: "/path/to/dir"

email
~~~~~

The value of ``email`` is needed if the pipeline is executed on a cluster,
which can use it to inform the person who started **uap** about status
changes of submitted jobs.

.. code-block:: yaml

    email: "your.name@mail.de"


tools
~~~~~

The ``tools`` block describes all programs needed during the execution of the
*uap**.

.. code-block:: yaml

    tools:
        # you don't have to specify a path if the tool can be found in $PATH
        cat:
            path: cat 
            version: "--version"
        # you have to specify a path if the tool can not be found in $PATH
        some-tool:
            path: /path/to/some-tool
            version: "--version"

steps
~~~~~

The ``steps`` block is the core of the analysis file, because it defines the
order in which the different steps of the analysis are executed.
Each step must have a unique name.
Therefore you should give each step a descriptive name followed by
a blank and the step type enclosed in parentheses.

There are two different types of steps:

1. **source steps** are used to enter data into the analysis, meaning they have no
   predecessor step they depend on.
2. **processing steps** depend upon one or more predecessor steps and create some 
   output that can be used by successor steps.
   
All available steps are described in detail in the steps documentation: 
:doc:`steps`.

Example configurations for various source steps are shown below:

.. code-block:: yaml

    # sources steps
    steps:
        # fastq_source provides a number of fastq.gz files as pipeline input
        casava_output (fastq_source):
            # a glob pattern
            pattern: /home/kaempf/Projects/RNAseq_Jurkats+BaP/data/
            group: (Sample_COPD_\d+)_R[12]-head.fastq.gz
            indices: indices.csv
            paired_end: yes

        # run_folder_sources
        fc1 (run_folder_source):
            path: /data/bioinf/projects/data/Jurkats_BaP_Transcriptome/130108_SN928_0083_AD11VNACXX_Keep/
            paired_end: yes
        fc2 (run_folder_source):
            path: /data/bioinf/projects/data/Jurkats_BaP_Transcriptome/130108_SN928_0084_BC0UT2ACXX_Keep/
            paired_end: yes
            
        # raw_file_source can provide any filesystem file as pipeline input
        mapped_reads (raw_file_source):
            path: data/H3K4me3_GCCAAT_L001_001.dup_rm.sam.gz
            sha1: 835779504aa63f80c9e1008f93f554269d0ec506
            
        # raw_url_source can provide any downloadable file as pipeline input
        gencode (raw_url_source):
            url: ftp://ftp.sanger.ac.uk/pub/gencode/release_15/gencode.v15.annotation.gtf.gz
            sha1: 9b272fde8bca544e6cd8621ddeec55aa09cf7a05



How-To Create Your First **uap** Analysis
=========================================

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


Example configurations
----------------------

There is currently no step implemented to execute Illuminas CASAVA pipeline, which 
converts BCL files to FASTQ files. Therefore all example configurations
begin with a source step that relies on the availability of fastq.gz files.

General sequencing analysis steps
********************************* 

Every analysis of high-throughput sequencing results starts with some basic
steps. Irrespective of sequencing RNA or DNA, given a reference genome
exists.

1. Get the sequencing reads as input (most likely fastq.gz)
2. Remove adapter sequences from your sequencing reads
3. Align the sequencing reads onto the refernce genome

After these steps are finished a lot of different analysis could be applied on
the data. Furtheron example configurations for often used analyses are shown.
The enumeration of steps show continues as if the basic steps were already
performed.


RNAseq analysis
***************


Differential expression
~~~~~~~~~~~~~~~~~~~~~~~

RNAseq analysis often aims at the discovery of differentially expressed
(known) transcripts. Therefore mappped reads for at least two different samples
have to be available.

4. Get annotation set (for e.g. genes, transcripts, ...)
5. Count the number of reads overlapping the annotation
6. Perform statistical analysis, based on counts 

Assemble novel transcripts
~~~~~~~~~~~~~~~~~~~~~~~~~~

As the publicly available annotations, e.g. from GENCODE, are probably not
complete, the assembly of novel transcripts from RNAseq data is another task one
would perform to invetsigate the transcriptome.


ChIPseq analysis
****************

ChIPseq analysis aims at the discovery of genomic loci at which protein(s) of
interest were bound. The experiment is an enrichment procedure using specific
antibodies. The enrichment detection is normally performed by so called peak
calling programs.

4. Get negative control
5. Peak calling


Prepare UCSC genome browser tracks
**********************************

The conversion of sequencing data into an format that can be displayed by the
UCSC genome browser is needed in almost all sequencing projects.


Interacting with a pipeline
===========================

Once the project is set up, there are several scripts which can be used to 
execute and monitor the pipeline. 
All scripts have a couple of properties in common:

* On startup, the configuration is read, tools are checked, input files are 
  collected, and all tasks are calculated. 
  If any of these steps fails, the script will print an error message with 
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
carried out and submits the jobs to a GridEngine cluster by calling ``qsub``. 
Dependencies are passed to ``qsub`` via the ``-hold_jid`` option, which means 
that jobs that depend on other jobs won't get scheduled until their 
dependencies have been satisfied. 
The file ``qsub-template.sh`` is used to submit jobs, with ``#{ }`` fields 
being substituted with appropriate values.

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

Post-mortem pipeline analysis
=============================
    
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


Extending rnaseq-pipeline
=========================


Implement your own steps
------------------------

The provided pipeline can be easily extended by implementing new steps and
sources. Therefore one does need some basic python programming skills. To add a
new processing step, a single Python file must be placed in ``include/step``
which defines a class with a constructor and two functions. The constructor
(``__init__``) checks for the availability of required tools and tells the
pipeline which connections this step expects (``in/``) and which it provides
(``out/``). The first of the functions  (``setup_runs``) is used for planning all
jobs based on a list of input files or runs and possibly additional information
from previous steps and the second function (``execute``) is used to execute a
specific job. The basic scaffold is shown below.

.. code-block:: python

    import sys
    from abstract_step import *
    import pipeline
    import re
    import process_pool
    import yaml
    
    class Macs14(AbstractStep):
        
        # the constructor
        def __init__(self, pipeline):
            super(Macs14, self).__init__(pipeline)

            # define in and out connections the strings have to start with 'in/'
            # or 'out/'
            self.add_connection('in/something')
            self.add_connection('out/tag1')
            self.add_connection('out/tag2')
            ...
    
            self.require_tool('cat4m')
            self.require_tool('pigz')
            ...

        # all checks of options and input values should be done here
        def setup_runs(self, complete_input_run_info, connection_info):
            # a hash containing information about this step
            output_run_info = {}

            # analyze the complete_input_run_info hash provided by the pipeline
            for step_name, step_input_info in complete_input_run_info.items():
                for input_run_id, input_run_info in step_input_info.items():
                   # assemble your output_run_info
                   # output_run_info has to look like this
                   output_run_info:
                       run_id_1:
                           "output_files":
                               tag1:
                                   output_file_1: [input_file_1, input_file_2, ...]
                                   output_file_2: [input_file_1, input_file_2, ...]
                               tag2:
                                   output_file_3: [input_file_1, input_file_2, ...]
                                   output_file_4: [input_file_1, input_file_2, ...]
                           "info":
                               ...
                           more:
                               ...
                           keys:
                               ...
                       run_id_2:
                           ...

            return output_run_info
        
        # called to actually launch the job (run_info is the hash returned from
        # setup_runs)
        def execute(self, run_id, run_info):
    
            with process_pool.ProcessPool(self) as pool:
                with pool.Pipeline(pool) as pipeline:
                    # assemble the steps pipline here
                    pipeline.append(...)
                    ...
                    # finally launch it
                    pool.launch(...)

The code shown above is the framework for a new step. The most essential part is
the hash returned by setup_runs(), here called ``output_run_info``.

:``run_id``:
    It has to be the unique name of a run (obviously, because its a key value).
    ``output_run_info`` can contain multiple ``run_id`` hashes.

:``"output_files"``:
    This is the only hash key that has to have a fix name. This is used to link
    input to output files.

:``tag[12]``:
    Every ``tag`` has to match ``\w+$`` in the string ``'out/tag'``, which was
    given to ``self.add_connection('out/tag')``. This can be any string, but it
    has to match with the last part of the connection string.

:``output_file_\d``:
    Each ``tag`` has to contain at least one such key. It has to be the name of
    the output file produced by the connection ``'out/tag'``. The value of this
    has to be a list of related input files. The list can have any number of
    entries even zero. Multiple ``output_file_\d`` can rely on the same set of
    input files.

Also very important is to understand the concept of *connections*. They provide
input files prior steps created already. The names of the connections can be
arbitrarily chosen, but should **not** describe the file format but more general
terms. For example an ``out/alignment`` can provide gzipped SAM or BAM files. So
you have to check in setup runs for the file type provided by a connection and
react accordingly. Inspect ``complete_input_run_info`` to find out what your
step gets as input.

Best practices
**************

There are a couple of things which should be kept in mind when implementing new 
steps or modifying existing steps:

* Make sure errors already show up in ``setup_runs`` instead of ``execute``.
  Therefore look out for things that may fail in ``setup_runs``. Stick to *fail
  early, fail often*. That way errors show up before submitting jobs to the
  cluster and wasting precious cluster waiting time is avoided. 
* Use the ``info`` entry in the returned ``output_run_info`` structure to pass
  information gathered in ``setup_runs`` to ``execute``.
* Likewise, make sure that the tools you'll need in ``execute`` are available.
  Check for the availability of tools within the constructor ``__init__``.

.. code-block:: python
  
    # make sure tools are available
    self.require_tool('pigz')
    self.require_tool('cutadapt')
    
* Make sure your disk access is as cluster-friendly as possible (which 
  primarily means using large block sizes and preferably no seek operations). 
  If possible, use ``unix_pipeline`` to wrap your commands in ``pigz``, ``dd``,
  or ``cat4m`` with a large block size like 4 MB. 
  Although this is not possible in every case (for example when seeking 
  in files is involved), it is straightforward with tools that read a 
  continuous stream from ``stdin`` and write a continuous stream to 
  ``stdout``.



Add the new step to your configuration
--------------------------------------

To insert a new step in a pipeline it has to be added into the ``config.yaml``.



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
