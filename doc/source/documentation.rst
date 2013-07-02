..
  This is the documentation for rnaseq-pipeline. Please keep lines under
  80 characters if you can and start each sentence on a new line as it 
  decreases maintenance and makes diffs more readable.

.. title:: Remaining original documentation


Setup
=====

The repository can be obtained like this::

    $ git clone spechtm@bioinf1:/home/spechtm/rnaseq-pipeline.git

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

All tools which are used in the pipeline must be specified in the 
configuration file.
The pipeline determines and records their versions for future reference.

By default, version determination is simply attempted by calling the program
without command-line arguments.

If a certain argument is required, specify it in ``get_version``. 
If the tools does not return with an exit code of 0, find out which code it
is by typing ``echo $?`` into Bash and specify the exit code in ``exit_code``.
            

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
