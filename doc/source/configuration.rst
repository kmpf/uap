..
  This is the documentation for rnaseq-pipeline. Please keep lines under
  80 characters if you can and start each sentence on a new line as it 
  decreases maintenance and makes diffs more readable.

.. title:: Configuring **uap**

..
  This document aims to describe how to configure **uap**.

Configure a Analysis for **uap**
================================

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

