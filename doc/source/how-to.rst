..
  This is the documentation for uap. Please keep lines under 80 characters if
  you can and start each sentence on a new line as it decreases maintenance
  and makes diffs more readable.

.. title:: How-To First Analysis

..
  This document describes how to set-up your first **uap** analysis.

How-To Create Your First **uap** Analysis
=========================================

Next, edit ``config.sample.yaml`` and save it as ``config.yaml``. 
Although writing the configuration may seem a bit complicated, the trouble 
pays off later because further interaction with the pipeline is quite simple. 
Here is a sample configuration:

.. code-block:: yaml

    # This is the rnaseq-pipeline configuration file.

    destination_path: "/home/michael/test-pipeline/out"

    steps:
        fastq_source:
            pattern: /home/michael/test-pipeline/fastq/*.fastq.gz
            group: (Sample_COPD_\d+)_R[12]-head.fastq.gz
            indices: copd-barcodes.csv
            paired_end: yes
            
        cutadapt:
            _depends: fastq_source
            adapter-R1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC((INDEX))ATCTCGTATGCCGTCTTCTGCTTG
            adapter-R2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
            
    tools:
        cutadapt:
            path: /home/michael/Desktop/rnaseq-pipeline/tools/cutadapt-1.2.1/bin/cutadapt
            get_version: '--version'
        pigz:
            path: pigz
            get_version: '--version'
        head:
            path: head
            get_version: '--version'
        cat4m:
            path: ./tools/cat4m


Example configurations
----------------------

There is currently no step implemented to execute Illuminas CASAVA pipeline, which 
converts BCL files to FASTQ files. Therefore all example configurations
begin with a source step that relies on the availability of fastq.gz files.

General sequencing analysis steps
********************************* 

Every analysis of high-throughput sequencing results starts with some basic
steps. Irrespective of sequencing RNA or DNA, given a reference genome
exists.

1. Get the sequencing reads as input (most likely fastq.gz)
2. Remove adapter sequences from your sequencing reads
3. Align the sequencing reads onto the refernce genome

After these steps are finished a lot of different analysis could be applied on
the data. Furtheron example configurations for often used analyses are shown.
The enumeration of steps show continues as if the basic steps were already
performed.


RNAseq analysis
***************


Differential expression
~~~~~~~~~~~~~~~~~~~~~~~

RNAseq analysis often aims at the discovery of differentially expressed
(known) transcripts. Therefore mappped reads for at least two different samples
have to be available.

4. Get annotation set (for e.g. genes, transcripts, ...)
5. Count the number of reads overlapping the annotation
6. Perform statistical analysis, based on counts 

Assemble novel transcripts
~~~~~~~~~~~~~~~~~~~~~~~~~~

As the publicly available annotations, e.g. from GENCODE, are probably not
complete, the assembly of novel transcripts from RNAseq data is another task one
would perform to invetsigate the transcriptome.


ChIPseq analysis
****************

ChIPseq analysis aims at the discovery of genomic loci at which protein(s) of
interest were bound. The experiment is an enrichment procedure using specific
antibodies. The enrichment detection is normally performed by so called peak
calling programs.

4. Get negative control
5. Peak calling


Prepare UCSC genome browser tracks
**********************************

The conversion of sequencing data into an format that can be displayed by the
UCSC genome browser is needed in almost all sequencing projects.

