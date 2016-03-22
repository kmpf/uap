..
  This is the documentation for uap. Please keep lines under
  80 characters if you can and start each sentence on a new line as it 
  decreases maintenance and makes diffs more readable.
  
.. title:: uap -- Universal Analysis Pipeline

.. _uap--index:
#########################################################
uap -- Robust, Consistent, and Reproducible Data Analysis
#########################################################

****

  **uap** executes, controls and keeps track of the analysis of large data sets.
  It enables users to perform robust, consistent, and reprodcuible data analysis.
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
  But, it does also support the cluster engines |uge_link|/OGE/SGE and
  |slurm_link|.


*****************
Table of contents
*****************

.. toctree::
   :maxdepth: 3

   how-to
   installation
   configuration
   interaction
   extension
   software-design
   annotation
   steps
   troubleshooting
   api


*******
Remarks
*******

This documentation has been created using |sphinx_link|
and |rest_link|.

******************
Indices and tables
******************

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
