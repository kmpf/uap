..
  This is the documentation for rnaseq-pipeline. Please keep lines under
  80 characters if you can and start each sentence on a new line as it 
  decreases maintenance and makes diffs more readable.

.. title:: Configuration of uap

..
  This document aims to describe how to configure **uap**.

Configuration of **uap**
========================

**uap** is made to control the execution of data analyses which are defined
in `YAML <http://www.yaml.org/>`_ files.
Each file describes a complete analysis.
Further on these files are called analysis file(s).

The configuration files consist of four sections (let's just call them sections,
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

Sections of a Configuration File
********************************

Destination_path Section
~~~~~~~~~~~~~~~~~~~~~~~~

The value of ``destination_path`` is the directory where **uap** is going
to store the created files. It is possible to use a different directory for
volatile files (see ).

.. code-block:: yaml

    destination_path: "/path/to/dir"

Email Section
~~~~~~~~~~~~~

The value of ``email`` is needed if the analysis is executed on a cluster,
which can use it to inform the person who started **uap** about status
changes of submitted jobs.

.. code-block:: yaml

    email: "your.name@mail.de"


Steps Section
~~~~~~~~~~~~~

The ``steps`` section is the core of the analysis file, because it defines when
steps are executed and how they depend on each other.
This section (technically it is a dictionary) contains a key for every step,
therefore each step must have a unique name.
There are two ways to name a step to allow multiple steps of the same type and
still ensure unique naming:

.. code-block:: yaml

    steps:
        # here, the step name is unchanged, it's a cutadapt step which is also
        # called 'cutadapt'
        cutadapt:
            ... # options following
            
        # here, we also insert a cutadapt step, but we give it a different name:
        # 'clip_adapters'
        clip_adapters (cutadapt):
            ... # options following
            
There are two different types of steps:

**Source Steps**
  They are used to provide files for the analysis. They have no dependencies
  and are usually the first steps of an analysis.

**Processing Steps**
  They depend upon one or more predecessor steps. Dependencies are defined via
  the ``_depends`` key which may either be ``null``, a step name, or a list of
  step names.

.. code-block:: yaml

    steps:
        # the source step which depends on nothing
        fastq_source:
            # ...
            
        # the first processing step, which depends on the source step
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

   
All available steps are described in detail in the steps documentation: 
:doc:`steps`.

.. _ToolsSection:
Tools Section
~~~~~~~~~~~~~

The ``tools`` section must list all programs needed during the execution of an
**uap** analysis.
**uap** determines and records their versions for future reference.

By default, version determination is simply attempted by calling the program
without command-line arguments.

If a certain argument is required, specify it in ``get_version``. 
If a tool does not exit with exit code 0, find out which code it is by typing
``echo $?`` into Bash and specify the exit code in ``exit_code``.



.. code-block:: yaml

    tools:
        # you don't have to specify a path if the tool can be found in $PATH
        cat:
            path: cat 
            get_version: "--version"
        # you have to specify a path if the tool can not be found in $PATH
        some-tool:
            path: /path/to/some-tool
            get_version: "--version"


Example Configurations
**********************

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

