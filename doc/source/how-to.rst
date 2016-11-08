..
  This is the documentation for uap. Please keep lines under 80 characters if
  you can and start each sentence on a new line as it decreases maintenance
  and makes diffs more readable.

.. title:: How-Tos

..
  This document describes how to set-up your first **uap** analysis.

.. _how-to:
##################
How-To Use **uap**
##################

At first, you need to install uap (see :ref:`installation_of_uap`).

***************************
Try Existing Configurations
***************************

After you have done that you need a working configuration file.
Example configurations are included in **uap**'s installation directory.
They are stored inside the ``example-configurations`` folder inside the
**uap** installation.
Go there and try::

  $ uap index_mycoplasma_genitalium_ASM2732v1_genome.yaml status

  Waiting tasks
  -------------
  [w] bowtie2_index/Mycoplasma_genitalium_index-download
  [w] bwa_index/Mycoplasma_genitalium_index-download
  [w] fasta_index/download
  [w] segemehl_index/Mycoplasma_genitalium_genome-download

  Ready tasks
  -----------
  [r] M_genitalium_genome/download

  tasks: 5 total, 4 waiting, 1 ready


Start your first **uap** analysis showcasing the controlled indexing of a
genome (arguably a tiny one)::

  $ uap index_mycoplasma_genitalium_ASM2732v1_genome.yaml status
  [uap] Set log level to ERROR
  [uap][ERROR]: index_mycoplasma_genitalium_ASM2732v1_genome.yaml: Destination path does not exist: genomes/bacteria/Mycoplasma_genitalium/

Oops, the :ref:`config_file_destination_path` does not exist.
Create it and start again::

  $ mkdir genomes/bacteria/Mycoplasma_genitalium/
  $ uap index_mycoplasma_genitalium_ASM2732v1_genome.yaml status

  Waiting tasks
  -------------
  [w] bowtie2_index/Mycoplasma_genitalium_index-download
  [w] bwa_index/Mycoplasma_genitalium_index-download
  [w] fasta_index/download
  [w] segemehl_index/Mycoplasma_genitalium_genome-download

  Ready tasks
  -----------
  [r] M_genitalium_genome/download

  tasks: 5 total, 4 waiting, 1 ready

If you still do get errors, please check that the tools defined in
``index_mycoplasma_genitalium_ASM2732v1_genome.yaml`` are available in your
environment (see :ref:`uap_config_tools_section`).

Go on and try some more example configurations (let's for now assume that all
tools are installed and configured correctly).
You want to create indexes of the human genome (hg19)::

  $ uap index_homo_sapiens_hg19_genome.yaml status
  [uap] Set log level to ERROR
  [uap][ERROR]: Output directory (genomes/animalia/chordata/mammalia/primates/homo_sapiens/hg19/chromosome_sizes) does not exist. Please create it.
  $ mkdir genomes/animalia/chordata/mammalia/primates/homo_sapiens/hg19/chromosome_sizes
  $ uap index_homo_sapiens_hg19_genome.yaml run-locally
  <Analysis starts>

*****************************
Create Your Own Configuration
*****************************

Although writing the configuration may seem a bit complicated, the trouble
pays off later because further interaction with the pipeline is quite simple.
The structure and content of the configuration files is very detailed described
on another page (see :ref:`configuration_of_uap`).
Here is a simple configuration:

.. code-block:: yaml

  Insert YAML here!

General Structure of Sequencing Analysis
========================================

Every analysis of high-throughput sequencing data evolves around some basic
tasks.
Irrespective of sequencing RNA or DNA.

1. Get the sequencing reads as input (most likely fastq.gz)
2. Remove adapter sequences from your sequencing reads
3. Align the sequencing reads onto the reference genome

The
After these steps are finished a lot of different analysis could be applied on
the data. Furtheron example configurations for often used analyses are shown.
The enumeration of steps show continues as if the basic steps were already
performed.


RNAseq analysis
---------------


Differential expression
^^^^^^^^^^^^^^^^^^^^^^^

RNAseq analysis often aims at the discovery of differentially expressed
(known) transcripts. Therefore mappped reads for at least two different samples
have to be available.

4. Get annotation set (for e.g. genes, transcripts, ...)
5. Count the number of reads overlapping the annotation
6. Perform statistical analysis, based on counts

Assemble novel transcripts
^^^^^^^^^^^^^^^^^^^^^^^^^^

As the publicly available annotations, e.g. from GENCODE, are probably not
complete, the assembly of novel transcripts from RNAseq data is another task one
would perform to invetsigate the transcriptome.


ChIPseq analysis
----------------

ChIPseq analysis aims at the discovery of genomic loci at which protein(s) of
interest were bound. The experiment is an enrichment procedure using specific
antibodies. The enrichment detection is normally performed by so called peak
calling programs.

4. Get negative control
5. Peak calling


Prepare UCSC genome browser tracks
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The conversion of sequencing data into an format that can be displayed by the
UCSC genome browser is needed in almost all sequencing projects.

