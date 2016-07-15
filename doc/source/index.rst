..
  This is the documentation for uap. Please keep lines under
  80 characters if you can and start each sentence on a new line as it 
  decreases maintenance and makes diffs more readable.
  
.. title:: uap -- Universal Analysis Pipeline

.. _uap--index:

#########################################################
uap -- Robust, Consistent, and Reproducible Data Analysis
#########################################################

  **uap** executes, controls and keeps track of the analysis of large data sets.
  It enables users to perform robust, consistent, and reproducible data analysis.
  **uap** encapsulates the usage of (bioinformatic) tools and handles data flow
  and processing during an analysis.
  Users can use predefined or self-made analysis steps to create custom analysis.
  Analysis steps encapsulate best practice usages for bioinformatic software
  tools.
  Although **uap**  focuses on the analysis of high-throughput sequencing data
  it can be extended to handle any analysis.
  But its plugin architecture allows users to add functionality.
  This would enable any kind of large data analysis.

**Usage:**

  **uap** is a command-line tool, implemented in Python, and runs under
  GNU/Linux.
  It takes a user-defined configuration file, which describes the analysis, as
  input.
  **uap** interacts with the analysis via subcommands.

**Supported Platforms:**

  **uap** runs natively on Unix-like operating systems.
  But, it does also support the cluster engines 
  |uge_link|, |oge_link| and |slurm_link|.

*********************
Important Information
*********************

The **uap** installation **does not** include all necessary tools for the data
analysis.
It expects that the required tools are **already installed**.

The recommended workflow to analyse data with **uap** is:

1. Install **uap** (see :doc:`installation`)
2. Optionally: Extend **uap** by adding new steps (see :doc:`extension`)
3. Write a configuration file to setup the analysis (see
   :doc:`configuration`)
4. Start the analysis locally (see :ref:`run-locally <uap-run-locally>`) or
   submit it to the cluster (see
   :ref:`submit-to-cluster <uap-submit-to-cluster>`)
5. Follow the progress of the analysis (see :ref:`status <uap-status>`)
6. Optionally: share your extensions with others (send a pull request via github)

A finished analysis leaves the user with:

* *The original input files* (which are, of course, left untouched).
* *The experiment-specific configuration file*
  (see :doc:`configuration`).
  You should keep this configuration file for later reference and you could
  even make it publicly available along with your input files for anybody to
  re-run the entire data analysis or parts thereof.
* *The result files and comprehensive annotations of the analysis*
  (see :doc:`annotation`).
  These files are located at the destination path defined in the configuration
  file.

************
Core aspects
************

Robust Data Analysis:
=====================

* Data is processed in temporary location.
  If and only if ALL involved processes exit gracefully, the output files are
  copied to the final output directory.
* Processing can be aborted and continued from the command line at any time.
  Failures during data processing do not lead to unstable state of analysis.
* The final output directory names are suffixed with a hashtag which is based
  on the commands executed to generate the output data.
  Data is not easily overwritten and this helps to check for necessary
  recomputations.
* Fail fast, fail often.
  During initiation **uap** controls the availability of all required tools,
  calculates the paths of all files of the analysis, and controls if these files
  exist.
  If any problems are encountered the user is requested to fix them.

Consistency:
============

* Steps and files are defined in a directed acyclic graph (DAG).
  The DAG defines dependencies between in- and output files.
* Prior to any execution the dependencies between files are calculated.
  If a file is newer or an option for a calculation has changed all dependent
  files are marked for recalculation.

Reproducibility:
================

* Comprehensive annotations are written to the output directories.
  They allow for later investigation of errors or review of executed commands.
  They contain also versions of used tool, required runtime, memory and CPU
  usage, etc.

Usability:
==========

* Single configuration file describes entire data analysis work-flow.
* Command-line tool to interacts with uap. It can be used to execute, monitor, and analyse defined work-flows.

*****************
Table of contents
*****************

.. toctree::
   :maxdepth: 3

   installation
   how-to
   configuration
   interaction
   extension
   software-design
   annotation
   steps
   api


*******
Remarks
*******

This documentation has been created using |sphinx_link| and |rest_link|.

******************
Indices and tables
******************

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. |uge_link| raw:: html
 
   <a href="http://www.univa.com/products/" target="_blank">UGE</a>

.. |oge_link| raw:: html

   <a href="http://www.univa.com/oracle" target="_blank">OGE/SGE</a>

.. |slurm_link| raw:: html
      
   <a href="http://slurm.schedmd.com/" target="_blank">SLURM</a>

.. |sphinx_link| raw:: html
 
    <a href="http://www.sphinx-doc.org/" target="_blank">sphinx</a>

.. |rest_link| raw:: html
 
    <a href="http://docutils.sourceforge.net/rst.html" target="_blank">reStructuredText</a>
