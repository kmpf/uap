..
  This is the documentation for rnaseq-pipeline. Please keep lines under
  80 characters if you can and start each sentence on a new line as it 
  decreases maintenance and makes diffs more readable.
  
.. title:: rnaseq-pipeline


Introduction
============

The **rnaseq-pipeline** package is a framework to configure, run, and control
high-throughput sequencing analyses.

The aim of this data processing pipeline is to enable robust and straightforward
bioinformatics data evaluation. It is implemented in Python, runs under
GNU/Linux and can be controlled from the command-line interface. Although the
primary focus is the evaluation of RNASeq data, its design allows for a variety
of other applications.

**Table of contents**

.. toctree::
   :maxdepth: 2

   documentation
   steps
   api


General usage
-------------

This package *does not* provide a number of tools which are downloaded and
installed system-wide to provide certain functioniality.
The intention of this system is to provide a robust and traceable framework
for data evaluation in scientific experiments which uses other tools and
manages individual data processing steps and their inter-dependencies.
    
The recommended workflow for running a data evaluation for an experiment is 
as follows:

1. Check-out the rnaseq-pipeline repository via Git.
2. Setup the project by writing the configuration file.
3. Add steps or other functionality as needed (optional).
4. Run the pipeline.
5. Have your changes (if there are any) merged back into the main repository,
   to the advantage of the scientific community.

This leaves you with:

* Your original input files, which are left untouched.
* The experiment-specific pipeline repository.  
  You should keep this repository for later reference and you could even
  make it publicly available along with your input files for anybody to
  re-run the entire data evaluation or parts thereof.
* The output directory containing all output files and comprehensive 
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
* Interaction with the pipeline happens through a handful of scripts which 
  are used to monitor the state of the pipeline and execute individual or all 
  remaining steps.

Design
------

The central part of the pipeline is its definition of the steps which are to 
be carried out.
Steps are organized in a dependency graph (a directed acyclic graph) -- every 
step may have one or more parent steps, which may in turn have other parent 
steps, and so on.
Steps without parents are usually sources which provide source files, for
example FASTQ files with the raw sequences obtained from the sequencer,
genome sequence databases or annotation tracks.

Each step defines a number of runs and each run represents a piece of the
entire data evaluation, typically at the level of a single sample.
A certain *run* of a certain *step* is called a *task*.
While the steps only describe what needs to be done on a very abstract level,
it is through the individual runs of each step that a pipeline-wide list of 
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

