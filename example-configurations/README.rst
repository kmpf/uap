Example workflows
==============

We provide 7 example configurations for workflows that can either be used directly with **uap** or
modified for your own purpose.

**Questions??**
- Why do we need the following two configuration file *index_mycoplasma_genitalium_ASM2732v1_genome.yaml*, it is very similar to 01.b

**ToDos!!**
- remove those lines from the worfklow configurations that are not required (commented addtional samples, files, etc.)
- rename the workflows to more intuitive files


01 Download input data and genome information
====

Before you can use one of the workflows, you need to download the example data.
Fortunately, you can do this also by using uap:

01.a ChIP-Seq data
----

This configuration downloads the example data for the ChIP-Seq data analysis workflow (see section
02, below) and saves the gzipped .fastq files into *example-out/2007-Barski_et_al_download*
directory.

**Configuration file:** 2007-CD4+_T_Cell_ChIPseq-Barski_et_al_download.yaml

**Steps:**

- download fastq files via ftp

01.b Genome data and mapping index files
----

This configuration file downloads chromosome 21 of the human genome (hg19) and creates the index
files for the tools used in the ChIP-Seq data analysis workflow (see section 02, below).

**Configuration file:** index_homo_sapiens_hg19_genome.yaml

**Steps:**

- download sequence of human chr21 (hg19)
- create index files for 
  - bowtie2
  - BWA
  -  fasta indices

01.c RNA CaptureSeq data
----

This configuration downloads the sequencing data analyzed in: *Targeted sequencing for gene
discovery and quantification using RNA CaptureSeq*, Mercer et al., Nature protocols, 2014. The data
is stored in the directory *example-out/2014-Mercer_et_al_download* and the workflow of the
respective data analysis can be found in section 03 below.

**Configuration file:** 2014-RNA_CaptureSeq-Mercer_et_al_download.yaml

**Steps:**

- download fastq and supplementary files via ftp and http
- unzip files
- compare secure hashes

01.d Gencode annotation
----

This configuration downloads the GENCODE annotation.

**Configuration file:** download_human_gencode_release.yaml

**Steps:**
- download .gtf annotation of hg19 genes via ftp
- download lncRNA .gtf annotation via ftp

02  ChIP-Seq data analysis
====

This workflow implements a High-Resolution Profiling of Histone Methylations in the Human Genome of
the data from *Barski et al., Cell 2007*. 

**Configuration file:** 2007-CD4+_T_Cell_ChIPseq-Barski_et_al.yaml

**Note:** 

If you want to run this workflow directly, it is required that you downloaded all needed data
previously. You can use the workflows 01.a and 01.b for this purpose.
If you want to modify this workflow for your own purpose, it might not be necessary to download the
example data.

**Steps:**

- read input data (fastq files)
- merge fastq files for each sample
- quality control (*fastqc* and *fastx quality stats*)
- adapter trimming (*cutadapt*) + QC
- read mapping onto genome (*bowtie*, *bwa*, *TopHat2*)
- sorting of alignments (*samtools*)
- mark duplicates (*picard tools*)
- peak calling (*MACS2*)

03 RNA CaptureSeq data analysis
====

This configuration repeats the analysis published in: *Targeted sequencing for gene discovery and
quantification using RNA CaptureSeq*, Mercer et al., Nature Protocols, 2014.

**Configuration file:** 2014-RNA_CaptureSeq-Mercer_et_al.yaml

**Note:** 

If you want to run this workflow directly, it is required that you download all needed data
previously. You can use the workflows 01.c for downloading the data and you need to modify workflow
01.b to download not only chromosome 21 but all chromosomes of the human genome.
In addition you need to download the GENCODE annotation file. You can use workflow 01.d for that.
If you want to modify this workflow for your own purpose, it might not be necessary to download the
example data.

**Steps:**
- read input data (fastq files)
- quality control (*fastqc* and *fastx quality stats*) 
- read mapping onto genome (*TopHat2*)
- sorting of alignments (*samtools*)
- assemble new transcripts (*cufflinks*)
- count reads overlapping genes (*htseq-count*)
