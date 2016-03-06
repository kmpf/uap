..
  This is the documentation for uap. Please keep lines under
  80 characters if you can and start each sentence on a new line as it 
  decreases maintenance and makes diffs more readable.
  
.. title:: uap -- Universal Analysis Pipeline


uap -- Robust, Consistent, and Reproducible Data Analysis
=========================================================

**uap** executes, controls and keeps track of the analysis of large data sets.
It enables users to perform robust, consistent, and reprodcuible data analysis.
**uap** is a command-line tool, implemented in Python, and runs under
GNU/Linux.
A single configuration file is required for an entire analysis.

Its main focus is the analysis of high-throughput sequencing data.
But, its plugin architecture allows users to add functionality and adapt it for
any kind of large data analysis.

General Information
-------------------

The **uap** installation *does not* include all necessary tools for the data
analysis.
It expects that the required tools are *already installed*.
**uap** encapsulates the usage of (bioinformatic) tools and handles data flow
and processing during an analysis.

The recommended workflow to analyse data with **uap** is:

1. Install **uap** (see :ref:`installation_of_uap`)
2. Optionally: Extend **uap** by adding new steps (see :ref:`extending_uap`)
3. Write a configuration file to setup the analysis (see
   :ref:`configuration_of_uap`)
4. Start the analysis locally (see :ref:`uap_run_locally`) or submit it to the
   cluster (see :ref:`uap_submit_to_cluster`)
5. Follow the progress of the analysis (see :ref:`uap_status`)
6. Share your extensions with others (send a pull request via github)

When the analysis is finished, you are left with:

* *The original input files*, which are, of course, left untouched
* *The experiment-specific configuration file* (see :ref:`configuration_of_uap`)
  You should keep this configuration file for later reference and you could even
  make it publicly available along with your input files for anybody to re-run
  the entire data analysis or parts thereof.
* *The destination path containing the output files and comprehensive 
  annotation files of the analysis* (see :ref:`annotations`).
  The annotation files contain detailed information about every output file.
  Also, the Git SHA1 hash of the **uap** repository at the time of
  data processing is included.
  The executed commands are listed.
  Annotation contains information about inter-process streams and output files,
  including SHA1 checksums, file sizes, and line counts as well.

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
  
**Reproduceability:**

* Comprehensive annotations are written to the output directories, allowing 
  for later investigation about what exactly happened.
      
**Simplicity:**

* The entire processing pipeline is described via a configuration file. 
  Steps are defined in a directed acyclic graph (DAG).
* Interaction with the pipeline happens through a single command-line tool which 
  can be used to execute and monitor the analysis.

Design
------

**uap** is designed as a plugin architecture.
The plugins are  where **uap** itself controls
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


Table of contents
=================

.. toctree::
   :maxdepth: 2

   how-to
   installation
   configuration
   interaction
   extension
   steps
   api



Remarks
=======

This documentation has been created using `Sphinx <http://sphinx-doc.org/>`_
and `reStructuredText <http://docutils.sourceforge.net/rst.html>`_.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

