..
  This is the documentation for uap. Please keep lines under
  80 characters if you can and start each sentence on a new line as it 
  decreases maintenance and makes diffs more readable.
  
.. title:: uap -- Universal Analysis Pipeline

.. _uap--index:

#########################################################
uap -- Robust, Consistent, and Reproducible Data Analysis
#########################################################

**Authors:**

Christoph Kämpf, Michael Specht, Sven-Holger Puppel, Alexander Scholz,
Gero Doose, Kristin Reiche, Jana Hertel, Jörg Hackermüller

**Description:**

  **uap** executes, controls and keeps track of the analysis of large data sets.
  It enables users to perform robust, consistent, and reproducible data analysis.
  **uap** encapsulates the usage of (bioinformatic) tools and handles data flow
  and processing during an analysis.
  Users can use predefined or self-made analysis steps to create custom analysis.
  Analysis steps encapsulate best practice usages for bioinformatic software
  tools.
  **uap** focuses on the analysis of high-throughput sequencing (HTS) data.
  But its plugin architecture allows users to add functionality, such that
  it can be used for any kind of large data analysis.

**Usage:**

  **uap** is a command-line tool, implemented in Python.
  It requires a user-defined configuration file, which describes the analysis,
  as input.
  
**Supported Platforms:**

  * Unix-like operating systems.
  * High Performance Compute (HPC) cluster systems such as |uge_link|,
    |oge_link| and |slurm_link|.
  * see :doc:`platforms` for detailed information
    
**Important Information**

  **uap** does **NOT** include all tools necessary for the data analysis.
  It expects that the required tools are **already installed**.

*****************
Table of contents
*****************

.. toctree::
   :maxdepth: 3

   introduction
   how-to
   installation
   configuration
   interaction
   extension
   software-design
   platforms
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
