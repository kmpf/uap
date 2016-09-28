
Steps
=====



Source steps
------------

Source steps provide input files for the pipeline, such as RNA sequences.

Run folder source
~~~~~~~~~~~~~~~~~

Here's an example:

.. code-block:: yaml

    - run_folder_source: { path: in }

This source looks for fastq.gz files in
``[path]/Unaligned/Project_*/Sample_*`` and pulls additional information from
CSV sample sheets it finds.
It also makes sure that index information for all samples is coherent and
unambiguous.

FASTQ source
~~~~~~~~~~~~

Here's an example:

.. code-block:: yaml

    - fastq_source:
        pattern: /data/original-fastq/&#42;.fastq.gz
        group: (Sample_COPD_\d+)_R[12].fastq.gz
        indices: copd-barcodes.csv

Input files are collected as defined by ``pattern`` and grouped into samples
according to ``group``, which is a regular expression.
All groups defined in the regex ``(  )`` are used to construct the sample
name, here it is used to declare that both R1 and R2 files belong to the
same sample.
Indices are read from the CSV file specified by ``indices``.

..
    .. automodule:: abstract_source

    .. autoclass:: AbstractSource
        :members:

*(detailed step descriptions to follow...)*

..
    Miscellaneous
    -------------

    Head
    ~~~~

    .. autosimpleclass:: head.Head

    Preprocessing
    -------------

    Adapter clipping
    ~~~~~~~~~~~~~~~~

    Cutadapt
    ^^^^^^^^

    .. autosimpleclass:: cutadapt.Cutadapt

    Fix cutadapt
    ^^^^^^^^^^^^

    .. autosimpleclass:: fix_cutadapt.FixCutadapt

    Aligners
    --------

    Segemehl
    ~~~~~~~~

    .. autosimpleclass:: segemehl.Segemehl
