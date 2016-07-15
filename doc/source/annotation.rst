..
  This is the documentation for uap. Please keep lines under 80 characters if
  you can and start each sentence on a new line as it decreases maintenance
  and makes diffs more readable.

.. title:: Results: Annotation Files

..
  This document aims to describe how to use **uap** via the command-line.

.. _annotation_files:

################
Annotation Files
################

The annotation files contain detailed information about every output file.
Also, the Git SHA1 hash of the **uap** repository at the time of
data processing is included.
The executed commands are listed.
Annotation contains information about inter-process streams and output files,
including SHA1 checksums, file sizes, and line counts as well.


Upon successful completion of a task, an extensive YAML-formatted annotation 
is placed next to the output files in a file called 
``.[task_id]-annotation.yaml``.
Also, for every output file, a symbolic link to this file is created:
``.[output_filename].annotation.yaml``.

Finally, the annotation is rendered via GraphViz, if available.
Rendering can also be done at a later time using annotations as input. (see :ref:`uap-render`)
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

known_paths
-----------

Contains information about all directories/files used during processing a run.
uap calculates the SHA1 hexdigest for each known file with the designation 'output' aka.
output/result files. 
