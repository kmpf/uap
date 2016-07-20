..
  This is the documentation for uap. Please keep lines under 80 characters if
  you can and start each sentence on a new line as it decreases maintenance
  and makes diffs more readable.

.. title:: How-Tos

..
  This document describes how to set-up your first **uap** analysis.

.. _how-to:

###################
Quick Start **uap**
###################

At first, you need to install **uap** (see :doc:`installation`).
After successfully finishing the installation of **uap** example data and
analysis can be found in the folder ``example-configurations``.

Let's jump head first into **uap** and have a look at some examples::

.. code-block:: bash

  $ cd <uap-path>/example-configurations/
  $ ls *.yaml
  2007-CD4+_T_Cell_ChIPseq-Barski_et_al_download.yaml
  2007-CD4+_T_Cell_ChIPseq-Barski_et_al.yaml
  ChIPseq-example-volatile.yaml
  ChIPseq-example.yaml
  ChIP_SEQ.yaml
  download_human_gencode_release.yaml
  index_homo_sapiens_hg19_genome.yaml
  index_mycoplasma_genitalium_ASM2732v1_genome.yaml
  RNAseq-example.yaml
  RNA_SEQ.yaml

..  $ uap index_mycoplasma_genitalium_ASM2732v1_genome.yaml status

These example configurations require different amounts of computational
resources.
Some example configurations download or work on small datasets and are
thus feasible for machines with limited resources.
Others require a very powerful stand-alone machine or a cluster system.
The examples are marked accordingly in the examples below.


Handle Genomic Data
-------------------

A usual analysis of High-Throughput Sequencing (HTS) data relies on different
publicly available data.
Most important is probably the genomic sequence of the species under
investigation.
That sequence is required to construct the indices (data structures used by 
read aligners).
Other publicly available data sets (such as reference annotations or the
chromosome sizes) migh also be required by the analysis.
The following configurations showcase how to get/generate that data:

*index_mycoplasma_genitalium_ASM2732v1_genome.yaml*
    Downloads the *Mycoplasma genitalium* genome, generates the indices for
    |bowtie2_link|, |bwa_link|, |segemehl_link|, and |samtools_link|.
    This workflow is quite fast because it uses the very small genome of
    *Mycoplasma genitalium*.

*index_homo_sapiens_hg19_genome.yaml*
    Downloads the *Homo sapiens* genome, generates the indices for
    |bowtie2_link|, |bwa_link|, |segemehl_link|, and |samtools_link|.
    This workflow requires substantial computational resources due to the
    size of the human genome.

*download_human_gencode_release.yaml*
    Downloads the human Gencode main annotation v24 and a subset for long
    non-coding RNA genes.
    This workflow only downloads files from the internetand and thus should
    work on any machine.

Let's have a look at the *Mycoplasma genitalium* example workflow by checking
its :ref:`uap_status`::

.. code-block:: bash

  $ cd <uap-path>/example-configurations/
  $ uap index_mycoplasma_genitalium_ASM2732v1_genome.yaml status
  [uap] Set log level to ERROR
  [uap][ERROR]: index_mycoplasma_genitalium_ASM2732v1_genome.yaml: Destination path does not exist: genomes/bacteria/Mycoplasma_genitalium/
  
Oops, the ``destination_path`` does not exist (see :ref:`config-file-destination-path`).
Create it and start again::

.. code-block:: bash

  $ mkdir -p genomes/bacteria/Mycoplasma_genitalium/
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

A list with all runs and their respective state should be displayed.
If the command still fails, please check that the tools defined in
``index_mycoplasma_genitalium_ASM2732v1_genome.yaml`` are available in your
environment (see :ref:`uap_config_tools_section`).
If you really want to download and index the genome tell **uap** to start
the workflow::

.. code-block:: bash

   $ uap index_mycoplasma_genitalium_ASM2732v1_genome.yaml run-locally

*2007-CD4+_T_Cell_ChIPseq-Barski_et_al_download.yaml*
    Download data published with the paper Barski et al., Cell (2007)

*2007-CD4+_T_Cell_ChIPseq-Barski_et_al.yaml*
    Perform some 
ChIPseq-example-volatile.yaml
ChIPseq-example.yaml
ChIP_SEQ.yaml
RNAseq-example.yaml
RNA_SEQ.yaml




The ``[w]`` stands for a waiting status of a task and the ``[r]`` stands for a ready status of a task. (see :doc:`interaction`)

Go on and try some more example configurations (let's for now assume that all
tools are installed and configured correctly).
You want to create indexes of the human genome (hg19)::

.. code-block:: bash

  $ uap index_homo_sapiens_hg19_genome.yaml status
  [uap] Set log level to ERROR
  [uap][ERROR]: Output directory (genomes/animalia/chordata/mammalia/primates/homo_sapiens/hg19/chromosome_sizes) does not exist. Please create it.
  $ mkdir genomes/animalia/chordata/mammalia/primates/homo_sapiens/hg19/chromosome_sizes
  $ uap index_homo_sapiens_hg19_genome.yaml run-locally
  <Analysis starts>

You get the idea.

*****************************
Create Your Own Configuration
*****************************

Although writing the configuration may seem a bit complicated, the trouble 
pays off later because further interaction with the pipeline is quite simple.
The structure and content of the configuration files is very detailed described
on another page (see :ref:`configuration-of-uap`).
Here is a simple configuration:

.. code-block:: yaml

  Insert YAML here!

General Structure of Sequencing Analysis
========================================

Every analysis of high-throughput sequencing (HTS) data evolves around some
basic tasks.
Irrespective of sequencing RNA or DNA.

1. Get the sequencing reads as input (most likely fastq.gz)
2. Remove adapter sequences from your sequencing reads
3. Align the sequencing reads onto the reference genome

After these steps are finished a lot of different analysis could be applied on
the data. Furtheron example configurations for often used analyses are shown.
The enumeration of steps continues as if the basic steps were already
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
would perform to investigate the transcriptome.


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


.. |bowtie2_link| raw:: html
      
   <a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml" target="_blank">bowtie2</a>

.. |bwa_link| raw:: html
      
   <a href="http://bio-bwa.sourceforge.net/" target="_blank">bwa</a>

.. |samtools_link| raw:: html
      
   <a href="http://www.htslib.org/" target="_blank">samtools</a>

.. |segemehl_link| raw:: html
      
   <a href="http://www.bioinf.uni-leipzig.de/Software/segemehl/" target="_blank">segemehl</a>

