..
  This is the documentation for uap. Please keep lines under 80 characters if
  you can and start each sentence on a new line as it decreases maintenance
  and makes diffs more readable.

.. title:: How-Tos

..
  This document describes how to set-up your first **uap** analysis.

.. _how-to:

*******************
Quick Start **uap**
*******************

At first, you need to install **uap** (see :doc:`installation`).
After successfully finishing the installation of **uap** example data and
analysis can be found in the folder ``example-configurations``.

Let's jump head first into **uap** and have a look at some examples::

  $ cd <uap-path>/example-configurations/
  $ ls *.yaml
  2007-CD4+_T_Cell_ChIPseq-Barski_et_al_download.yaml
  2007-CD4+_T_Cell_ChIPseq-Barski_et_al.yaml
  ChIPseq-example.yaml
  download_human_gencode_release.yaml
  index_homo_sapiens_hg19_genome.yaml
  index_mycoplasma_genitalium_ASM2732v1_genome.yaml
  RNAseq-example.yaml


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

``index_mycoplasma_genitalium_ASM2732v1_genome.yaml``
    Downloads the *Mycoplasma genitalium* genome, generates the indices for
    |bowtie2_link|, |bwa_link|, |segemehl_link|, and |samtools_link|.
    This workflow is quite fast because it uses the very small genome of
    *Mycoplasma genitalium*.

``index_homo_sapiens_hg19_genome.yaml``
    Downloads the *Homo sapiens* genome, generates the indices for
    |bowtie2_link|, |bwa_link|, and |samtools_link|.
    This workflow requires substantial computational resources due to the
    size of the human genome.
    The |segemehl_link| index creation is commented out due to its high
    memory consumption (~50-60 GB).
    Please make sure to only run it on well equipped machines.

``download_human_gencode_release.yaml``
    Downloads the human Gencode main annotation v24 and a subset for long
    non-coding RNA genes.
    This workflow only downloads files from the internetand and thus should
    work on any machine.

Let's have a look at the *Mycoplasma genitalium* example workflow by checking
its :ref:`uap_status`::

  $ cd <uap-path>/example-configurations/
  $ uap index_mycoplasma_genitalium_ASM2732v1_genome.yaml status
  [uap] Set log level to ERROR
  [uap][ERROR]: index_mycoplasma_genitalium_ASM2732v1_genome.yaml: Destination path does not exist: genomes/bacteria/Mycoplasma_genitalium/
  
Oops, the ``destination_path`` does not exist (see :ref:`config-file-destination-path`).
Create it and start again::

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
A run is always in one of these states:

* ``[r]eady``                                                                                                                                                    
* ``[w]aiting``                                                                                                                                                  
* ``[q]ueued``                                                                                                                                                   
* ``[e]xecuting``                                                                                                                                                
* ``[f]inished``

If the command still fails, please check that the tools defined in
``index_mycoplasma_genitalium_ASM2732v1_genome.yaml`` are available in your
environment (see :ref:`uap_config_tools_section`).
If you really want to download and index the genome tell **uap** to start
the workflow::

   $ uap index_mycoplasma_genitalium_ASM2732v1_genome.yaml run-locally

**uap** should have created a symbolic link named
``index_mycoplasma_genitalium_ASM2732v1_genome.yaml-out`` pointing to the 
``destination_path``.
The content should look something like that::

    $ tree --charset=ascii
    .
    |-- bowtie2_index
    |   |-- Mycoplasma_genitalium_index-download-cMQPtBxs
    |   |   |-- Mycoplasma_genitalium_index-download.1.bt2
    |   |   |-- Mycoplasma_genitalium_index-download.2.bt2
    |   |   |-- Mycoplasma_genitalium_index-download.3.bt2
    |   |   |-- Mycoplasma_genitalium_index-download.4.bt2
    |   |   |-- Mycoplasma_genitalium_index-download.rev.1.bt2
    |   |   `-- Mycoplasma_genitalium_index-download.rev.2.bt2
    |   `-- Mycoplasma_genitalium_index-download-ZsvbSjtK
    |       |-- Mycoplasma_genitalium_index-download.1.bt2
    |       |-- Mycoplasma_genitalium_index-download.2.bt2
    |       |-- Mycoplasma_genitalium_index-download.3.bt2
    |       |-- Mycoplasma_genitalium_index-download.4.bt2
    |       |-- Mycoplasma_genitalium_index-download.rev.1.bt2
    |       `-- Mycoplasma_genitalium_index-download.rev.2.bt2
    |-- bwa_index
    |   `-- Mycoplasma_genitalium_index-download-XRyj5AnJ
    |       |-- Mycoplasma_genitalium_index-download.amb
    |       |-- Mycoplasma_genitalium_index-download.ann
    |       |-- Mycoplasma_genitalium_index-download.bwt
    |       |-- Mycoplasma_genitalium_index-download.pac
    |       `-- Mycoplasma_genitalium_index-download.sa
    |-- fasta_index
    |   `-- download-HA439DGO
    |       `-- Mycoplasma_genitalium.ASM2732v1.fa.fai
    |-- M_genitalium_genome
    |   `-- download-5dych7Xj
    |-- Mycoplasma_genitalium.ASM2732v1.fa
    |-- segemehl_index
    |   |-- Mycoplasma_genitalium_genome-download-2UKxxupJ
    |   |   |-- download-segemehl-generate-index-log.txt
    |   |   `-- Mycoplasma_genitalium_genome-download.idx
    |   `-- Mycoplasma_genitalium_genome-download-zgtEpQmV
    |       |-- download-segemehl-generate-index-log.txt
    |       `-- Mycoplasma_genitalium_genome-download.idx
    `-- temp

Congratulation you've finished your first **uap** workflow!

Go on and try to run some more workflows.
Most examples require the human genome so you might turn your head towards the
``index_homo_sapiens_hg19_genome.yaml`` workflow from her::

  $ uap index_homo_sapiens_hg19_genome.yaml status
  [uap] Set log level to ERROR
  [uap][ERROR]: Output directory (genomes/animalia/chordata/mammalia/primates/homo_sapiens/hg19/chromosome_sizes) does not exist. Please create it.
  $ mkdir -p genomes/animalia/chordata/mammalia/primates/homo_sapiens/hg19/chromosome_sizes
  $ uap index_homo_sapiens_hg19_genome.yaml run-locally
  <Analysis starts>

Again you need to create the output folder (you get the idea).


Sequencing Data Analysis
------------------------

Now that you possess the genome sequences, indices, and annotations let's have
a look at some example analysis.

General Steps
^^^^^^^^^^^^^

The analysis of high-throughput sequencing (HTS) data usually start with some
basic steps.

1. Conversion of the raw sequencing data to, most likely, fastq(.gz) files
2. Removal of adapter sequences from the sequencing reads
3. Alignment of the sequencing reads onto the reference genome

These basic steps can be followed up with a lot of different analysis steps.
The following analysis examples illustrate how to perform the basic as well as
some more specific steps.

RNAseq Example Workflow
^^^^^^^^^^^^^^^^^^^^^^^

``RNAseq-example.yaml``
    This workflow includes the aforementioned basic steps, followed by the
    use of |tophat2_link| to map split reads.
    After mapping |htseq_count_link| is used to count the reads mapping to
    every exon from the annotation.
    A |cufflinks_link| step for *de novo* transcript assembly is also included.
    The usage of **|segemehl_link|** is disabled by default.
    But it can be enabled and combined with |cufflinks_link| *de novo*
    transcript assembly employing our **s2c** python script.
    **This workflow is not going to work, because the initial data set is
    to small.**

Differential expression
"""""""""""""""""""""""

RNAseq analysis often aims at the discovery of differentially expressed
(known) transcripts. Therefore mappped reads for at least two different samples
have to be available.

4. Get annotation set (for e.g. genes, transcripts, ...)
5. Count the number of reads overlapping the annotation
6. Perform statistical analysis, based on counts 

Assemble novel transcripts
""""""""""""""""""""""""""

The publicly available annotations, e.g. from GENCODE, might be not
complete.
That's why the assembly of novel transcripts from RNAseq data is a common
task.

ChIPseq Example Workflow
^^^^^^^^^^^^^^^^^^^^^^^^

``ChIPseq-example.yaml``
    This workflow includes the aforementioned basic steps, followed by the
    use of |bowtie2_link| to map the sequencing reads.
    After the mapping duplicate reads are removed using |picard_link|.
    And finally |macs2_link| is used to infer enriched regions aka. peak
    calling.
    **This workflow is not going to work, because the initial data set is
    to small.**

Peak Calling
""""""""""""

ChIPseq analysis aims at the discovery of genomic loci at which protein(s) of
interest were bound. The experiment is an enrichment procedure using specific
antibodies. The enrichment detection is normally performed by so called peak
calling programs.

4. Duplicate removal
5. Peak calling

Real-Life Data Set Example
^^^^^^^^^^^^^^^^^^^^^^^^^^

An example **uap** workflow is included in the two configuration files:

* ``2007-CD4+_T_Cell_ChIPseq-Barski_et_al_download.yaml``
* ``2007-CD4+_T_Cell_ChIPseq-Barski_et_al.yaml``


``2007-CD4+_T_Cell_ChIPseq-Barski_et_al_download.yaml``
    Downloads the data published with the paper |Barski_link|.

``2007-CD4+_T_Cell_ChIPseq-Barski_et_al.yaml``

************************
Create Your Own Workflow
************************

You finished to check out the examples?
Go and try to create your own workflow
If you are fine with what you saw 
Although writing the configuration may seem a bit complicated, the trouble 
pays off later because further interaction with the pipeline is quite simple.
The structure and content of the configuration files is very detailed described
on another page (see :ref:`configuration-of-uap`).


.. |Barski_link| raw:: html
   <a href="http://www.ncbi.nlm.nih.gov/pubmed/17512414" target="_blank">Barski <i>et al.</i>, Cell (2007)</a>

.. |bowtie2_link| raw:: html
      
   <a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml" target="_blank">bowtie2</a>

.. |bwa_link| raw:: html
      
   <a href="http://bio-bwa.sourceforge.net/" target="_blank">bwa</a>

.. |cufflinks_link| raw:: html
   
   <a href="" target="_blank">cufflinks</a>

.. |htseq_count_link| raw:: html
      
   <a href="http://www-huber.embl.de/users/anders/HTSeq/doc/count.html" target="_blank">htseq-count</a>

   .. |picard_link| raw:: html
      
   <a href="http://broadinstitute.github.io/picard/" target="_blank">Picard</a>

.. |samtools_link| raw:: html
      
   <a href="http://www.htslib.org/" target="_blank">samtools</a>

.. |segemehl_link| raw:: html
      
   <a href="http://www.bioinf.uni-leipzig.de/Software/segemehl/" target="_blank">segemehl</a>

.. |tophat2_link| raw:: html
      
   <a href="https://ccb.jhu.edu/software/tophat/index.shtml" target="_blank">tophat2</a>
