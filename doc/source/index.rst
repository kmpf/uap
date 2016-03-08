..
  This is the documentation for uap. Please keep lines under
  80 characters if you can and start each sentence on a new line as it 
  decreases maintenance and makes diffs more readable.
  
.. title:: uap -- Universal Analysis Pipeline

.. _uap--index

uap -- Robust, Consistent, and Reproducible Data Analysis
=========================================================

**uap** executes, controls and keeps track of the analysis of large data sets.
It enables users to perform robust, consistent, and reprodcuible data analysis.
Users can either combine predefined analysis steps to create custom analysis or
they can extend **uap** with their own analysis steps.
Steps are best practice usages for the encapsulated commands.

**uap** is a command-line tool, implemented in Python, and runs under
GNU/Linux.
It takes a user-defined configuration file, which describes the analysis, as
input.
**uap** works on the analysis via subcommands.

**uap**'s  main focus is the analysis of high-throughput sequencing data.
But, as already mentioned, its plugin architecture allows users to add
functionality.
This would enable any kind of large data analysis.

Important Information
---------------------

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

* *The original input files* (which are, of course, left untouched)
* *The experiment-specific configuration file* (see :ref:`configuration_of_uap`)
  You should keep this configuration file for later reference and you could even
  make it publicly available along with your input files for anybody to re-run
  the entire data analysis or parts thereof.
* *The output files and comprehensive annotations of the analysis* (see :ref:`annotations`).
  These files are stored in the destination path defined in the configuration file.

Core aspects
------------

**Robustness:**

* All steps write their output files to a temporary location. 
  Only if a step has completed successfully, the output files are copied to 
  the correct output directory.
* The output directory names are suffixed with a hashtag which is based on the
  commands executed to generate the output data.
* Processing can be aborted and continued from the command line at any time. 
  This way, cluster failures are less critical because output files do not
  get compromised.
* Errors are caught as early as possible. Tools are checked for availability, 
  and the entire processing pipeline is calculated in advance before 
  jobs are being started or submitted to a cluster.

**Consistency:**

* Steps and files are defined in a directed acyclic graph (DAG).
  The DAG defines dependencies between in- and output files.
* Prior to any execution the dependencies between files are calculated.
  If a file is newer or an option for a calculation has changed all dependent
  files are marked for recalculation.

**Reproducibility:**

* Comprehensive annotations are written to the output directories.
  They allow for later investigation of errors or review of executed commands.
  They contain also versions of used tool, required runtime, memory and CPU
  usage, etc.

**Simplicity:**

* The entire processing pipeline is described via a configuration file. 
* Interaction with the pipeline happens through a single command-line tool which 
  can be used to execute and monitor the analysis.


Table of contents
=================

.. toctree::
   :maxdepth: 2

   how-to
   installation
   configuration
   interaction
   extension
   software-design
   annotation
   steps
   api
   post-mortem


Remarks
=======

This documentation has been created using `Sphinx <http://sphinx-doc.org/>`_
and `reStructuredText <http://docutils.sourceforge.net/rst.html>`_.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

