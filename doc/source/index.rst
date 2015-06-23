..
  This is the documentation for uap. Please keep lines under
  80 characters if you can and start each sentence on a new line as it 
  decreases maintenance and makes diffs more readable.
  
.. title:: uap -- Universal Analysis Pipeline


Introduction
============

The **uap** package is a framework to configure, run, and control
large data multi-step analyses.
Its main focus is on the analysis of high-throughput sequencing data.

The aim of this data processing pipeline is to enable robust and straightforward
bioinformatics data evaluation.
It is implemented in Python, runs under GNU/Linux and can be controlled from the
command-line interface.
Although the primary focus is the evaluation of sequencing data, its design
allows for a variety of other applications.

**Table of contents**

.. toctree::
   :maxdepth: 2

   installation
   configuration
   interaction
   extension
   how-to
   steps
   api


General Information
-------------------

This package *does not* provide a number of tools which are downloaded and
installed system-wide to provide certain functioniality.
The intention of this system is to provide a robust and traceable framework
for data evaluation in scientific experiments which uses other tools and
manages individual data processing steps and their inter-dependencies.
    
The recommended workflow for running a data evaluation for an experiment is 
as follows:

1. Check-out the **uap** repository from `GitHub <https://github.com/kmpf/uap>`_.
2. Make **uap** globally available (as described )
2. Setup the analysis by writing a configuration file.
3. Add steps or other functionality as needed (optional).
4. Start the execution of the analysis.
5. Monitor the execution with the command-line tools.
5. Share your added functionalities (if there are any) with the scientific 
   community by merging them back into the main repository.

This leaves you with:

* Your original input files, which are, of course, left untouched.
* The experiment-specific configuration file.
  You should keep this configuration file for later reference and you could even
  make it publicly available along with your input files for anybody to re-run
  the entire data analysis or parts thereof.
* An output directory containing all output files and comprehensive 
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
* Interaction with the pipeline happens through a single command-line tool which 
  can be used to execute and monitor the analysis.

Design
------

**uap** is designed as a plugin architecture, where **uap** itself controls
the ordered execution of the plugged in so called steps.
Steps are organized in a dependency graph (a directed acyclic graph) -- every 
step may have one or more parent steps, which may in turn have other parent 
steps, and so on.
Steps without parents are usually sources which provide source files, for
example FASTQ files with the raw sequences obtained from the sequencer,
genome sequence databases or annotation tracks.

Each step defines a number of runs and each run represents a piece of the
entire data analysis, typically at the level of a single sample.
A certain *run* of a certain *step* is called a *task*.
While the steps only describe what needs to be done on a very abstract level,
it is through the individual runs of each step that a **uap** wide list of 
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

            
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

