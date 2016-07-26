..
  This is the documentation for uap. Please keep lines under
  80 characters if you can and start each sentence on a new line as it 
  decreases maintenance and makes diffs more readable.
  
.. title:: Introducing **uap**

*******************
Introducing **uap**
*******************

.. _uap-core-aspects:

Core aspects
============

Robustness
----------

* Data is processed in temporary location.
  If and only if ALL involved processes exited graceful, the output files are
  copied to the final output directory.
* The final output directory names are suffixed with a hashtag which is based
  on the commands executed to generate the output data.
  Data is not easily overwritten and this helps to check for necessary
  recomputations.
* Processing can be aborted and continued from the command line at any time.
  Failures during data processing do not lead to unstable state of analysis.
* Errors are reported as early as possible, fail-fast.
  Tools are checked for availability, and the entire processing pipeline is
  calculated in advance before jobs are being started or submitted to a cluster.

.. _uap-consistency:

Consistency
-----------

* Steps and files are defined in a directed acyclic graph (DAG).
  The DAG defines dependencies between in- and output files.
* Prior to any execution the dependencies between files are calculated.
  If a file is newer or an option for a calculation has changed all dependent
  files are marked for recalculation.

Reproducibility
---------------

* Comprehensive annotations are written to the output directories.
  They allow for later investigation of errors or review of executed commands.
  They contain also versions of used tool, required runtime, memory and CPU
  usage, etc.

Usability
---------

* Single configuration file describdes entire processing pipeline.
* Single command-line tool interacts with the pipeline.
  It can be used to execute, monitor, and analyse the pipeline.

.. _uap-recommended-workflow:

Recommended Workflow
====================

The recommended workflow to analyse data with **uap** is:

1. Install **uap** (see :doc:`installation`)
2. Optionally: Extend **uap** by adding new steps (see :doc:`extension`)
3. Write a configuration file to setup your analysis (see
   :doc:`configuration`)
4. Start the analysis locally (see :ref:`run-locally <uap-run-locally>`) or
   submit it to a cluster (see
   :ref:`submit-to-cluster <uap-submit-to-cluster>`)
5. Follow the progress of the analysis (see :ref:`status <uap-status>`)
6. Share your extensions with the public (send us a pull request via github)

A **finished** analysis leaves the user with:

* *The original input files* (which are, of course, left untouched).
* *The experiment-specific configuration file*
  (see :doc:`configuration`).
  You should keep this configuration file for later reference and you could
  even make it publicly available along with your input files for anybody to
  re-run the entire data analysis or parts thereof.
* *The output files and comprehensive annotations of the analysis*
(see :doc:`annotation`).
  These files are stored in the destination path defined in the configuration
  file.


.. |uge_link| raw:: html

   <a href="http://www.univa.com/products/" target="_blank">UGE</a>

.. |slurm_link| raw:: html

   <a href="http://slurm.schedmd.com/" target="_blank">SLURM</a>

.. |sphinx_link| raw:: html

   <a href="http://sphinx-doc.org/" target="_blank">Sphinx</a>

.. |rest_link| raw:: html

   <a href="http://docutils.sourceforge.net/rst.html" target="_blank">`reStructuredText</a>
