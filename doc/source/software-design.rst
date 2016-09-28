..
  This is the documentation for uap. Please keep lines under 80 characters if
  you can and start each sentence on a new line as it decreases maintenance
  and makes diffs more readable.

.. title:: Software Design

..
  This document aims to describe how to use **uap** via the command-line.

.. _software_design:
###############
Software Design
###############

**uap** is designed as a plugin architecture.
The plugins are internally called steps.
Two different types of steps exist, the source and processing steps.
Source steps are used to include data from outside the destination path (see
:ref:`config_file_destination_path`) into the analysis.
Processing steps are blueprints
Each step corresponds to the blueprint of a single data processing. where **uap** itself controls
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
