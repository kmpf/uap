Example workflows
#################

We provide eleven example configurations for workflows that can either be used
directly with **uap** or modified for your own purpose.

These four configurations download genome information and create derived data
structures required for the data analysis workflows:

- download_human_gencode_release_v19.yaml
- index_homo_sapiens_hg19_chr21.yaml
- index_homo_sapiens_hg19_genome.yaml
- index_mycoplasma_genitalium_ASM2732v1_genome.yaml

These three configurations download sequencing data for the data analysis
workflows:
    
- ChIPseq-data-download-full.yaml
- ChIPseq-data-download-short.yaml
- RNAseq-data-download.yaml

These four configurations implement the data analysis:
    
- ChIPseq-workflow-full.yaml
- ChIPseq-workflow-short.yaml
- RNAseq-workflow-full.yaml
- RNAseq-workflow-short.yaml


1. Download genome information and create indices for the data analysis workflows
=================================================================================

Most analysis require some kind of genome sequence and/or annotation data.

1.1 download_human_gencode_release_v19.yaml
-------------------------------------------

This configuration downloads the GENCODE annotation.

**Configuration file:** download_human_gencode_release.yaml

**Steps:**

- download .gtf annotation for hg19 genes via ftp
- download lncRNA .gtf annotation for hg19 via ftp


1.2 index_homo_sapiens_hg19_chr21.yaml
--------------------------------------

This configuration file downloads the FASTA sequence of chromosome 21 of the
human genome (version hg19) and creates the index files for the tools used in
the data analysis workflows (see section 3, below).

**Configuration file:** index_homo_sapiens_hg19_chr21.yaml

**Steps:**

- download sequence of human chr21 (hg19)
- create index files for:
  - bowtie2/tophat2
  - BWA
  - fasta indices

1.3 index_homo_sapiens_hg19_genome.yaml
---------------------------------------

This configuration file downloads the FASTA sequence of all human chromosomes
(version hg19) and creates the index files for the tools used in the data
analysis workflows (see section 3, below).

**Configuration file:** index_homo_sapiens_hg19_genome.yaml

**Steps:**

- download sequences of human chr1-chr22, chrX and chrY (hg19)
- merge all chromosme sequences
- create index files for:
  - bowtie2/tophat2
  - BWA
  - fasta indices
- download chromosome sizes file

1.4 index_mycoplasma_genitalium_ASM2732v1_genome.yaml
-----------------------------------------------------

This configuration file downloads the FASTA sequence of the Mycoplasma
genitalium genome and creates the index files for different read mapping tools.
*Due to the small genome size the creation of the segemehl index can be tested
with this workflow.*

**Configuration file:** index_mycoplasma_genitalium_ASM2732v1_genome.yaml

**Steps:**

- download genome sequence of Mycoplasma genitalium
- create index files for:
  - bowtie2/tophat2
  - BWA
  - segemehl
  - fasta indices


2. Download sequencing data for the data analysis workflows
===========================================================

The ChIPseq and RNAseq analysis (see section 3) are based on published data
sets. The following configurations download the required sequencing data.

2.1 ChIPseq-data-download-full.yaml & ChIPseq-data-download-short.yaml
----------------------------------------------------------------------

This configuration downloads the example data for the ChIP-Seq data analysis
workflow (see section 02, below) and saves the gzipped .fastq files into
*example-out/2007-Barski_et_al_download* directory.

**Configuration file:** ChIPseq-data-download-full.yaml

**Steps:**
- download fastq files with all ChIPseq data sets via ftp

**Configuration file:** ChIPseq-data-download-short.yaml

**Steps:**
- download fastq files with H3K4me1 and H3K4me3 ChIPseq data sets via ftp


2.2 RNAseq-data-download.yaml
-----------------------------

This configuration downloads the sequencing data analyzed in: *Targeted
sequencing for gene discovery and quantification using RNA CaptureSeq*, Mercer
et al., Nature protocols, 2014. The data is stored in the directory
*example-out/2014-Mercer_et_al_download* and the workflow of the respective data
analysis can be found in section 3.2 below.

**Configuration file:** RNAseq-data-download.yaml

**Steps:**
- download fastq and supplementary files via ftp and http
- unzip files
- compare secure hashes


3. Data analysis workflows
==========================

Most analysis require some kind of genome sequence and/or annotation data.

3.1 ChIPseq-workflow-full.yaml &  ChIPseq-workflow-short.yaml
-------------------------------------------------------------

This workflow implements a High-Resolution Profiling of Histone Methylations in
the Human Genome of the data from *Barski et al., Cell 2007*.

**Configuration file:**
- ChIPseq-workflow-short.yaml
- ChIPseq-workflow-full.yaml

**Note:**

If you want to run this workflow directly, it is required that you downloaded
all needed data previously. You can use the workflows in section 2.1 for this
purpose. If you want to modify this workflow for your own purpose, it might not
be necessary to download the example data.

**Steps:**

- read input data (fastq files)
- merge fastq files for each sample
- quality control (*fastqc* and *fastx quality stats*)
- adapter trimming (*cutadapt*) + QC
- read mapping onto genome (*bowtie*, *bwa*, *TopHat2*)
- sorting of alignments (*samtools*)
- mark duplicates (*picard tools*)
- peak calling (*MACS2*)


3.2 RNAseq-workflow-short.yaml & RNAseq-workflow-full.yaml
----------------------------------------------------------

This configuration repeats the analysis published in: *Targeted sequencing for
gene discovery and quantification using RNA CaptureSeq*, Mercer et al., Nature
Protocols, 2014.

**Configuration file:**
- RNAseq-workflow-short.yaml
- RNAseq-workflow-full.yaml

**Note:** 

If you want to run this workflow directly, it is required that you download all
needed data previously. You can use the workflow in section 2.2 for downloading
the data and you need to modify workflow 01.b to download not only chromosome 21 but all
chromosomes of the human genome. In addition you need to download the GENCODE
annotation file. You can use workflow 01.d for that. If you want to modify this
workflow for your own purpose, it might not be necessary to download the example
data.

**Steps:**

- RNAseq-workflow-short.yaml
 - read input data (fastq files)
 - quality control (*fastqc* and *fastx quality stats*) 
 - read mapping onto genome (*TopHat2*)
 - sorting of alignments (*samtools*)
 - assemble new transcripts (*cufflinks*)
 - count reads mapped by tophat2 overlapping genes (*htseq-count*)

- RNAseq-workflow-full.yaml (contains all the steps in RNAseq-workflow-short.yaml plus)
 - read mapping onto genome (*segemehl*)
 - make segemehl output compatible with cufflinks (*s2c*)
 - sort alignments by position (*samtools*)
 - assemble new transcripts (*cufflinks*)
 - count reads mapped by segemehl overlapping genes (*htseq-count*)
