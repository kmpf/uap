<!--
  This is the documentation for uap. Please keep lines under 80 characters if
  you can and start each sentence on a new line as it decreases maintenance
  and makes diffs more readable.
-->

<!--
  This document describes how to set-up your first **uap** analysis.
-->

# Quick Start

At first, you need to install **uap** (see [Installation](./installation.md)).
After successfully finishing the installation of **uap** example
analysis can be found in the folder
`example-configurations`.

Let's jump head first into **uap** and have a look at some examples:

```bash
  $ cd <uap-path>/example-configurations/
  $ ls *.yaml
  2007-CD4+_T_Cell_ChIPseq-Barski_et_al_download.yaml
  2007-CD4+_T_Cell_ChIPseq-Barski_et_al.yaml
  2014-RNA_CaptureSeq-Mercer_et_al_download.yaml
  2014-RNA_CaptureSeq-Mercer_et_al.yaml
  download_human_gencode_release.yaml
  index_homo_sapiens_hg19_genome.yaml
  index_mycoplasma_genitalium_ASM2732v1_genome.yaml
```

These example configurations differ in their usage of computational
resources.
Some example configurations download or work on small datasets and are
thus feasible for machines with limited resources.
Most examples can be extended by uncommenting additional steps.
This might change their computational requirements in such a way that a
very powerful stand-alone machine or a cluster system is required.
The examples are marked accordingly in the sections below.

**NOTE**: Before **computing an example on a cluster**, you need to uncomment
the [`cluster` section](./configuration.md#section-cluster) and adapt
the settings as required.
Please check also if the
[Cluster Configuration File](./configuration.md#cluster-configuration-file)
fits your cluster system.

**NOTE**: The examples contain information where users can obtain
**required external/bioinformatics tools**.
If **uap** fails due to a missing tool, please check the
provided URLs for installation instructions.

## Handle Genomic Data

A usual analysis of High-Throughput Sequencing (HTS) data relies on different
publicly available data.
Most important is probably the genomic sequence of the species under
investigation.
That sequence is required to construct the indices (data structures used by 
read aligners).
Other publicly available data sets (such as reference annotations or the
chromosome sizes) might also be required for an analysis.
The following configurations showcase how to get or generate that data:

### `index_mycoplasma_genitalium_ASM2732v1_genome.yaml`

Downloads the *Mycoplasma genitalium* genome, generates the indices for
[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml),
[bwa](http://bio-bwa.sourceforge.net/),
[segemehl](http://www.bioinf.uni-leipzig.de/Software/segemehl/), and
[samtools](http://www.htslib.org/).
This workflow is quite fast because it uses the very small genome of
*Mycoplasma genitalium*.

* Max. memory: ~0,5 GB
* Disk usage: ~20 MB
* Run time: minutes 
* Required tools:
  * [bwa](http://bio-bwa.sourceforge.net/)
  * [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
  * [curl](https://curl.haxx.se/)
  * [pigz](http://zlib.net/pigz/)
  * [samtools](http://www.htslib.org/)
  * [segemehl](http://www.bioinf.uni-leipzig.de/Software/segemehl/)

### `index_homo_sapiens_hg19_chr21.yaml`

Downloads chromosome 21 of the *Homo sapiens* genome, generates the indices for
[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml),
[bwa](http://bio-bwa.sourceforge.net/), and
[samtools](http://www.htslib.org/).
This minimal version should work just fine.
Users can uncomment steps to download the complete genome.
This would substantially increase the required computational resources.
The [segemehl](http://www.bioinf.uni-leipzig.de/Software/segemehl/) index
creation is commented out due to its high memory consumption (~50-60 GB), if
working with the whole genome.

* Max. memory: ~2 GB
* Disk usage: ~240 MB
* Run time: several minutes
* Required tools:
  * [bwa](http://bio-bwa.sourceforge.net/)
  * [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
  * [curl](https://curl.haxx.se/)
  * [fetchChromSizes](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes)
  * [pigz](http://zlib.net/pigz/)
  * [samtools](http://www.htslib.org/)
  * [segemehl](http://www.bioinf.uni-leipzig.de/Software/segemehl/) (if uncommented)


### `index_homo_sapiens_hg19_genome.yaml`

Downloads *Homo sapiens* chromosome 21, generates the indices for
[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml),
[bwa](http://bio-bwa.sourceforge.net/), and
[samtools](http://www.htslib.org/) (and
[segemehl](http://www.bioinf.uni-leipzig.de/Software/segemehl/) if
you uncomment it).
This workflow requires substantial computational resources due to the
size of the human genome.
The [segemehl](http://www.bioinf.uni-leipzig.de/Software/segemehl/) index
creation is commented out due to its high memory consumption.
Please make sure to only run it on well equipped machines.

* Required tools:
  * [bwa](http://bio-bwa.sourceforge.net/)
  * [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
  * [curl](https://curl.haxx.se/)
  * [fetchChromSizes](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes)
  * [pigz](http://zlib.net/pigz/)
  * [samtools](http://www.htslib.org/)
  * [segemehl](http://www.bioinf.uni-leipzig.de/Software/segemehl/) (if uncommented)
    
### `download_human_gencode_release.yaml`

Downloads the human Gencode main annotation v19 and a subset for long
non-coding RNA genes.
This workflow only downloads files from the internet and and thus should
work on any machine.

* Max. memory: depends on your machine
* Disk usage: ~1,2 GB
* Run time: depends on your internet connection
* Required tools:
  * [bwa](http://bio-bwa.sourceforge.net/)
  * [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
  * [curl](https://curl.haxx.se/)
  * [fetchChromSizes](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes)
  * [pigz](http://zlib.net/pigz/)
  * [samtools](http://www.htslib.org/)

Let's have a look at the *Mycoplasma genitalium* example workflow by checking
its [status](./interaction.md#subcommand-status):

```bash
  $ cd <uap-path>/example-configurations/
  $ uap index_mycoplasma_genitalium_ASM2732v1_genome.yaml status
  [uap][ERROR]: index_mycoplasma_genitalium_ASM2732v1_genome.yaml: Destination path does not exist: genomes/bacteria/Mycoplasma_genitalium/
```

Oops, the `destination_path` does not exist (see
[section `destination_path`](./configuration.md#section-destination_path)).
Create it and start again:

```bash
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
```

A list with all runs and their respective state should be displayed.
A run is always in one of these states:

* ``[r]eady``
* ``[w]aiting``
* ``[q]ueued``
* ``[e]xecuting``
* ``[f]inished``
* ``[c]hanged``
* ``[b]ad``

If the command still fails, please check that the tools defined in
`index_mycoplasma_genitalium_ASM2732v1_genome.yaml` are available in your
environment (see
[section `tools`](./configuration.md#section-tools)).
If you really want to download and index the genome tell **uap** to start
the workflow:

```bash
   $ uap index_mycoplasma_genitalium_ASM2732v1_genome.yaml run-locally
```

**uap** should have created a symbolic link named
`index_mycoplasma_genitalium_ASM2732v1_genome.yaml-out` pointing to the 
`destination_path`.
The content should look something like that:

```bash
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
```

Congratulation you've finished your first **uap** workflow!

Go on and try to run some more workflows.
Most examples require the human genome so you might turn your head towards the
`index_homo_sapiens_hg19_genome.yaml` workflow from her:

```bash
  $ uap index_homo_sapiens_hg19_genome.yaml status
  [uap][ERROR]: Output directory (genomes/animalia/chordata/mammalia/primates/homo_sapiens/hg19/chromosome_sizes) does not exist. Please create it.
  $ mkdir -p genomes/animalia/chordata/mammalia/primates/homo_sapiens/hg19/chromosome_sizes
  $ uap index_homo_sapiens_hg19_genome.yaml run-locally
  <Analysis starts>
```

Again you need to create the output folder (you get the idea).
Be aware that by default only the smallest chromosome, chromsome 21, is
downloaded and indexed.
This reduces required memory and computation time.
You can uncomment the download steps for the other chromosomes and the index
for the complete genome will be created.

## Sequencing Data Analysis

Now that you possess the genome sequences, indices, and annotations let's have
a look at some example analysis.

### General Steps

The analysis of high-throughput sequencing (HTS) data usually start with some
basic steps.

1. Conversion of the raw sequencing data to, most likely, fastq(.gz) files
2. Removal of adapter sequences from the sequencing reads
3. Alignment of the sequencing reads onto the reference genome

These basic steps can be followed up with a lot of different analysis steps.
The following analysis examples illustrate how to perform the basic as well as
some more specific steps.

### RNAseq Example -- Reanalysing Data from [Mercer *et al.*, Nature Protoc. (2014)](https://www.ncbi.nlm.nih.gov/pubmed/24705597)

RNAseq analysis often aims at the discovery of differentially expressed
(known) transcripts. Therefore mappped reads for at least two different samples
have to be available.

A. Differential Expression Analysis

   4. Get annotation set (for e.g. genes, transcripts, ...)
   5. Count the number of reads overlapping the annotation
   6. Perform statistical analysis, based on counts 

Another common analysis performed with RNAseq data is the identification of
novel tarnscripts. This approach is useful to identify tissue-specific
transcipts.
      
B. *De novo* Transcript Assembly
   
   4. Apply transcript assembly tool on mapped reads

`2014-RNA_CaptureSeq-Mercer_et_al_download.yaml`
    Downloads the data published in the paper
    [Mercer *et al.*, Nature Protoc. (2014)](https://www.ncbi.nlm.nih.gov/pubmed/24705597).

    :Max. memory: ~? GB
    :Disk usage: ~12 GB
    :Run time: minutes (depending on your internet connection)

    Required tools:

        * [curl](https://curl.haxx.se/)
        * [pigz](http://zlib.net/pigz/)

`2014-RNA_CaptureSeq-Mercer_et_al.yaml`
    The downloaded FASTQ files get analysed by [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and
    [FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/).
    The reads are afterwards mapped to the human genome with [tophat2](https://ccb.jhu.edu/software/tophat/index.shtml).
    The mapped reads are afterwards sorted by position using [samtools](http://www.htslib.org/).
    [htseq-count](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html) is used to count the mapped reads for every exon of
    the annotation.
    [cufflinks](http://cufflinks.cbcb.umd.edu/) is used to perform *de novo* transcript assembly.
    The usage of [segemehl](http://www.bioinf.uni-leipzig.de/Software/segemehl/) is **disabled** by default.
    But it can be enabled and combined with [cufflinks](http://cufflinks.cbcb.umd.edu/) *de novo*
    transcript assembly employing our **s2c** python script.

    :Max. memory: ~? GB
    :Disk usage: ~3 GB
    :Run time: several hours

        * [cufflinks](http://cufflinks.cbcb.umd.edu/)
        * [cutadapt](https://github.com/marcelm/cutadapt)
        * [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
        * [FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)
        * [htseq-count](http://www-huber.embl.de/users/anders/HTSeq/doc/count.html)
        * [pigz](http://zlib.net/pigz/)
        * [samtools](http://www.htslib.org/)
        * [segemehl](http://www.bioinf.uni-leipzig.de/Software/segemehl/) (if uncommented)
        * [tophat2](https://ccb.jhu.edu/software/tophat/index.shtml)
               
.. NOTE:: Before computing ``2014-RNA_CaptureSeq-Mercer_et_al.yaml``
          please make sure that, the following examples were executed:

          - `index_homo_sapiens_hg19_genome.yaml`
          - `download_human_gencode_release.yaml`

### ChIPseq Example: Reanalysing Data from [Barski *et al.*, Cell (2007)](http://www.ncbi.nlm.nih.gov/pubmed/17512414)

ChIPseq analysis aims at the discovery of genomic loci at which protein(s) of
interest were bound. The experiment is an enrichment procedure using specific
antibodies. The enrichment detection is normally performed by so called peak
calling programs. The data is prone to duplicate reads from PCR due to relatively
low amounts of input DNA. So these steps follow the basic ones:

4. Duplicate removal
5. Peak calling

The analysis of data published in the paper
[Barski *et al.*, Cell (2007)](http://www.ncbi.nlm.nih.gov/pubmed/17512414) is
contained in these files:

.. _example_barski_download:

### `2007-CD4+_T_Cell_ChIPseq-Barski_et_al_download.yaml`

Downloads the data published in the paper
[Barski *et al.*, Cell (2007)](http://www.ncbi.nlm.nih.gov/pubmed/17512414).

Max. memory: ~? GB
Disk usage: ~17 GB
Run time: depends on your internet connection

Downloads the data published in the paper
[Barski *et al.*, Cell (2007)](http://www.ncbi.nlm.nih.gov/pubmed/17512414).

Required tools:

* [curl](https://curl.haxx.se/)
* [pigz](http://zlib.net/pigz/)

.. _example_barski:

### `2007-CD4+_T_Cell_ChIPseq-Barski_et_al.yaml`

At first the downloaded FASTQ files are grouped by sample.
All files per sample are merged. 
Sequencing quality is controlled by
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
and [FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/).
Adapter sequences are removed from the reads before they are mapped to 
the human genome.
Reads are mapped with
[bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml),
[bwa](http://bio-bwa.sourceforge.net/), and
[tophat2](https://ccb.jhu.edu/software/tophat/index.shtml).
Again mapping with
[segemehl](http://www.bioinf.uni-leipzig.de/Software/segemehl/)
is disabled by default due to its high resource requirements.
Library complexity is estimated using
[preseq](http://smithlabresearch.org/software/preseq/).
After the mapping duplicate reads are removed using
[Picard](http://broadinstitute.github.io/picard/).
Finally enriched regions are detected with
[MACS2](https://github.com/taoliu/MACS).

Max. memory: ~? GB
Disk usage: ~51 GB
Run time: ~several hours (on a cluster), ~1 day (on a single machine)

Required tools:

* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [bwa](http://bio-bwa.sourceforge.net/)
* [cutadapt](https://github.com/marcelm/cutadapt)
* [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)
* [MACS2](https://github.com/taoliu/MACS)
* [Picard](http://broadinstitute.github.io/picard/)
* [pigz](http://zlib.net/pigz/)
* [preseq](http://smithlabresearch.org/software/preseq/)
* [samtools](http://www.htslib.org/)
* [segemehl](http://www.bioinf.uni-leipzig.de/Software/segemehl/)
* [tophat2](https://ccb.jhu.edu/software/tophat/index.shtml)

HINT: The usage of [Picard](http://broadinstitute.github.io/picard/) can differ
a lot between systems.

On Ubuntu systems it can be called like this:

```bash
$ picard-tools --version
```

If you use it as recommended at
[Picard](http://broadinstitute.github.io/picard/), it is called like this:

```bash
$ java -jar /path/to/picard.jar -h
```

Please check how to use it on your system and adjust the example
configuration accordingly (see [`tools`](./configuration.md#section-tools)).
            
NOTE: Before computing `2007-CD4+_T_Cell_ChIPseq-Barski_et_al.yaml`
please make sure that, the following examples were executed:

- `index_homo_sapiens_hg19_genome.yaml`
- `download_human_gencode_release.yaml`

## Create Your Own Workflow

You finished to check out the examples?
Go and try to create your own workflow
If you are fine with what you saw.
Although writing the configuration may seem a bit complicated, the trouble 
pays off later because further interaction with the pipeline is quite simple.
The structure and content of the configuration files is very detailed described
on another page (see
[Analysis Configuration File](./configuration.md#analsis-configuration-file)).
