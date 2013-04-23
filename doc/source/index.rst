Introduction
============

The aim of this data processing pipeline is to enable simple and robust bioinformatics data evaluation.

**Simplicity:**

* The entire processing pipeline is described via a config file. Step are defined in a tree, and output files are written into a directory structure mirroring this tree.
* To add a new processing step, a single Python file must be placed in ``include/step`` which defines a class with two functions, one for planning all jobs based on a list of input files or runs and possibly additional information from previous steps and another function for running a specific job.

**Robustness:**

* All steps write their output files to a temporary location (a fact which a step is not aware of). Only if a step has completed successfully, the output files are copied to the correct output directory.
* The output directory names are suffixed with a eight-character hashtag which mirrors the options specified for the step.
* Processing can be aborted and continued from the command line at any time. This way, cluster failures are less critical.
* Comprehensive annotations are written to the output directories, allowing for later investigation.
* Errors are caught as early as possible. Tools are checked for availability, the entire processing pipeline is calculated in advance.

A pipeline is defined by two aspects:

* the steps it carries out, with dependencies defined via a tree
* its input samples

The combination of *steps* and *samples* result in a list of *tasks*, which can be executed sequentially or can be submitted to a cluster.

.. IMPORTANT:: The design decision that steps are defined as a tree instead of a full directed acyclic graph means that a step cannot have more than one direct parent, like a directory in a file system cannot have more than one parent directory. This means that a step cannot use the output of two different steps as its input.

Setup
=====

After cloning the repository, run the bootstrapping script to create the required Python environment (which will be located in ``./python_env/``)::

    $ ./bootstrap.sh

There's no harm in accidentally running this script multiple times.

The configuration file
----------------------

Next, edit ``config.sample.yaml`` and save it as ``config.yaml``. Although writing the configuration may seem a bit complicated, it pays off later because further interaction with the pipeline is quite simple. Here is a sample configuration::

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

* ``sources`` - there can be multiple sources of different types:

  * run folders
  * plain fastq.gz files with additional information
  
* ``destination_path`` - this is where result files, annotations and temporary files are written to
* ``steps`` - defines the processing step arranged in a tree
* ``tools`` - defines all tools used in the pipeline and how to determine their versions (for later reference)
* ``email`` - when submitting jobs on a cluster, messages will be sent to this email address (nobody@example.com by default)

The pipeline
------------

.. toctree::
   :maxdepth: 2
   
.. automodule:: pipeline

.. autoclass:: Pipeline
    :members:
    
Sources
~~~~~~~
    
.. automodule:: abstract_source

.. autoclass:: AbstractSource
    :members:
    
Steps
-----

.. automodule:: abstract_step
    :members:

.. autoclass:: AbstractStep
    :members:
    
Available steps
~~~~~~~~~~~~~~~
    
.. automodule:: head
.. autoclass:: Head
    :members:

.. automodule:: break
.. autoclass:: Break
    :members:

.. automodule:: cutadapt
.. autoclass:: Cutadapt
    :members:

.. automodule:: fix_cutadapt
.. autoclass:: FixCutadapt
    :members:

.. automodule:: segemehl
.. autoclass:: Segemehl
    :members:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

