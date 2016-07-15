###############
Available steps
###############

************
Source steps
************

.. index:: bcl2fastq_source

bcl2fastq_source
================

**Connections:**
  - Output Connection:
    
    - 'out/configureBcl2Fastq_log_stderr'
    - 'out/make_log_stderr'
    - 'out/sample_sheet'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      bcl2fastq_source [style=filled, fillcolor="#fce94f"];
      out_0 [label="configureBcl2Fastq_log_stderr"];
      bcl2fastq_source -> out_0;
      out_1 [label="make_log_stderr"];
      bcl2fastq_source -> out_1;
      out_2 [label="sample_sheet"];
      bcl2fastq_source -> out_2;
   }    

**Options:**
  - **adapter-sequence** (str, optional)
    
  - **adapter-stringency** (str, optional)
    
  - **fastq-cluster-count** (int, optional)
    
  - **filter-dir** (str, optional)
    
  - **flowcell-id** (str, optional)
    
  - **ignore-missing-bcl** (bool, optional)
    
  - **ignore-missing-control** (bool, optional)
    
  - **ignore-missing-stats** (bool, optional)
    
  - **input-dir** (str, required) -- file URL
    
  - **intensities-dir** (str, optional)
    
  - **mismatches** (int, optional)
    
  - **no-eamss** (str, optional)
    
  - **output-dir** (str, optional)
    
  - **positions-dir** (str, optional)
    
  - **positions-format** (str, optional)
    
  - **sample-sheet** (str, required)
    
  - **tiles** (str, optional)
    
  - **use-bases-mask** (str, optional) -- Conversion mask characters:- Y or y: use- N or n: discard- I or i: use for indexingIf not given, the mask will be guessed from theRunInfo.xml file in the run folder.For instance, in a 2x76 indexed paired end run, themask *Y76,I6n,y75n* means: "use all 76 bases from thefirst end, discard the last base of the indexing read,and use only the first 75 bases of the second end".
    
  - **with-failed-reads** (str, optional)
    
**Required tools:** configureBclToFastq.pl, make, mkdir, mv

This step provides input files which already exists and therefore creates no tasks in the pipeline.

.. index:: fastq_source

fastq_source
============

The FastqSource class acts as a source for FASTQ files. This source creates a
run for every sample.

Specify a file name pattern in *pattern* and define how sample names should
be determined from file names by specifyign a regular expression in *group*.

Sample index barcodes may specified by providing a filename to a CSV file
containing the columns *Sample_ID* and *Index* or directly by defining a
dictionary which maps indices to sample names.

**Connections:**
  - Output Connection:
    
    - 'out/first_read'
    - 'out/second_read'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      fastq_source [style=filled, fillcolor="#fce94f"];
      out_0 [label="first_read"];
      fastq_source -> out_0;
      out_1 [label="second_read"];
      fastq_source -> out_1;
   }    

**Options:**
  - **first_read** (str, required) -- Part of the file name that marks all files containing sequencing data of the first read. Example: 'R1.fastq' or '_1.fastq'
    
  - **group** (str, optional) -- A regular expression which is applied to found files, and which is used to determine the sample name from the file name. For example, ``(Sample_\d+)_R[12].fastq.gz``, when applied to a file called ``Sample_1_R1.fastq.gz``, would result in a sample name of ``Sample_1``. You can specify multiple capture groups in the regular expression.
    
  - **indices** (str/dict, optional) -- path to a CSV file or a dictionary of sample_id: barcode entries.
    
  - **paired_end** (bool, required) -- Specify whether the samples are paired end or not.
    
  - **pattern** (str, optional) -- A file name pattern, for example ``/home/test/fastq/Sample_*.fastq.gz``.
    
  - **sample_id_prefix** (str, optional) -- This optional prefix is prepended to every sample name.
    
  - **sample_to_files_map** (dict/str, optional) -- A listing of sample names and their associated files. This must be provided as a YAML dictionary.
    
  - **second_read** (str, required) -- Part of the file name that marks all files containing sequencing data of the second read. Example: 'R2.fastq' or '_2.fastq'
    
This step provides input files which already exists and therefore creates no tasks in the pipeline.

.. index:: fetch_chrom_sizes_source

fetch_chrom_sizes_source
========================

**Connections:**
  - Output Connection:
    
    - 'out/chromosome_sizes'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      fetch_chrom_sizes_source [style=filled, fillcolor="#fce94f"];
      out_0 [label="chromosome_sizes"];
      fetch_chrom_sizes_source -> out_0;
   }    

**Options:**
  - **path** (str, required) -- directory to move file to
    
  - **ucsc-database** (str, required) -- Name of UCSC database e.g. hg38, mm9
    
**Required tools:** cp, fetchChromSizes

This step provides input files which already exists and therefore creates no tasks in the pipeline.

.. index:: raw_file_source

raw_file_source
===============

**Connections:**
  - Output Connection:
    
    - 'out/raw'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      raw_file_source [style=filled, fillcolor="#fce94f"];
      out_0 [label="raw"];
      raw_file_source -> out_0;
   }    

**Options:**
  - **group** (str, optional) -- A regular expression which is applied to found files, and which is used to determine the sample name from the file name. For example, `(Sample_\d+)_R[12].fastq.gz``, when applied to a file called ``Sample_1_R1.fastq.gz``, would result in a sample name of ``Sample_1``. You can specify multiple capture groups in the regular expression.
    
  - **pattern** (str, optional) -- A file name pattern, for example ``/home/test/fastq/Sample_*.fastq.gz``.
    
  - **sample_id_prefix** (str, optional)
    
  - **sample_to_files_map** (dict/str, optional) -- A listing of sample names and their associated files. This must be provided as a YAML dictionary.
    
This step provides input files which already exists and therefore creates no tasks in the pipeline.

.. index:: raw_file_sources

raw_file_sources
================

The RawFileSources class acts as a temporary fix to get files into the pipeline.
This source creates a run for every sample.

Specify a file name pattern in *pattern* and define how sample names should be
determined from file names by specifyign a regular expression in *group*.

**Connections:**
  - Output Connection:
    
    - 'out/raws'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      raw_file_sources [style=filled, fillcolor="#fce94f"];
      out_0 [label="raws"];
      raw_file_sources -> out_0;
   }    

**Options:**
  - **group** (str, required) -- A regular expression which is applied to found files, and which is used to determine the sample name from the file name. For example, ``(Sample_\d+)_R[12].fastq.gz``, when applied to a file called ``Sample_1_R1.fastq.gz``, would result in a sample name of ``Sample_1``. You can specify multiple capture groups in the regular expression.
    
  - **paired_end** (bool, required) -- Specify whether the samples are paired end or not.
    
  - **pattern** (str, required) -- A file name pattern, for example ``/home/test/fastq/Sample_*.fastq.gz``.
    
  - **sample_id_prefix** (str, optional) -- This optional prefix is prepended to every sample name.
    
This step provides input files which already exists and therefore creates no tasks in the pipeline.

.. index:: raw_url_source

raw_url_source
==============

**Connections:**
  - Output Connection:
    
    - 'out/raw'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      raw_url_source [style=filled, fillcolor="#fce94f"];
      out_0 [label="raw"];
      raw_url_source -> out_0;
   }    

**Options:**
  - **filename** (str, optional) -- local file name of downloaded file
    
  - **hashing-algorithm** (str, optional) -- hashing algorithm to use
    
    - possible values: 'md5', 'sha1', 'sha224', 'sha256', 'sha384', 'sha512'
    
  - **path** (str, required) -- directory to move downloaded file to
    
  - **secure-hash** (str, optional) -- expected secure hash of downloaded file
    
  - **uncompress** (bool, optional) -- Shall the file be uncompressed after downloading
    
  - **url** (str, required) -- file URL
    
**Required tools:** compare_secure_hashes, cp, curl, dd, mkdir, pigz

This step provides input files which already exists and therefore creates no tasks in the pipeline.

.. index:: run_folder_source

run_folder_source
=================

This source looks for fastq.gz files in
``[path]/Unaligned/Project_*/Sample_*`` and pulls additional information
from CSV sample sheets it finds. It also makes sure that index information
for all samples is coherent and unambiguous.

**Connections:**
  - Output Connection:
    
    - 'out/first_read'
    - 'out/second_read'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      run_folder_source [style=filled, fillcolor="#fce94f"];
      out_0 [label="first_read"];
      run_folder_source -> out_0;
      out_1 [label="second_read"];
      run_folder_source -> out_1;
   }    

**Options:**
  - **first_read** (str, required) -- Part of the file name that marks all files containing sequencing data of the first read. Example: '_R1.fastq' or '_1.fastq'
    
  - **paired_end** (bool, required)
    
  - **path** (str, required)
    
  - **project** (str, required)
    
  - **second_read** (str, required) -- Part of the file name that marks all files containing sequencing data of the second read. Example: 'R2.fastq' or '_2.fastq'
    
This step provides input files which already exists and therefore creates no tasks in the pipeline.

****************
Processing steps
****************

.. index:: bam_to_genome_browser

bam_to_genome_browser
=====================

**Connections:**
  - Input Connection:
    
    - 'in/alignments'
  - Output Connection:
    
    - 'out/alignments'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      bam_to_genome_browser [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> bam_to_genome_browser;
      out_1 [label="alignments"];
      bam_to_genome_browser -> out_1;
   }    

**Options:**
  - **bedtools-bamtobed-color** (str, optional)
    
  - **bedtools-bamtobed-tag** (str, optional)
    
  - **bedtools-genomecov-3** (bool, optional)
    
  - **bedtools-genomecov-5** (bool, optional)
    
  - **bedtools-genomecov-max** (int, optional)
    
  - **bedtools-genomecov-report-zero-coverage** (bool, required)
    
  - **bedtools-genomecov-scale** (float, optional)
    
  - **bedtools-genomecov-split** (bool, required)
    
  - **bedtools-genomecov-strand** (str, optional)
    
    - possible values: '+', '-'
    
  - **chromosome-sizes** (str, required)
    
  - **output-format** (str, required)
    
    - possible values: 'bed', 'bigBed', 'bedGraph', 'bigWig'
    
  - **trackline** (dict, optional)
    
  - **trackopts** (dict, optional)
    
**Required tools:** bedGraphToBigWig, bedToBigBed, bedtools, dd, mkfifo, pigz

**CPU Cores:** 8

.. index:: bowtie2

bowtie2
=======

Bowtie2 is an ultrafast and memory-efficient tool for aligning sequencing
reads to long reference sequences. It is particularly good at aligning reads
of about 50 up to 100s or 1,000s of characters, and particularly good at
aligning to relatively long (e.g. mammalian) genomes. Bowtie 2 indexes the
genome with an FM Index to keep its memory footprint small: for the human
genome, its memory footprint is typically around 3.2 GB. Bowtie 2 supports
gapped, local, and paired-end alignment modes.

|bowtie2_link|

typical command line:

.. code-block:: bash

    bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r>} -S [<hit>]

**Connections:**
  - Input Connection:
    
    - 'in/first_read'
    - 'in/second_read'
  - Output Connection:
    
    - 'out/alignments'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      bowtie2 [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> bowtie2;
      in_1 [label="second_read"];
      in_1 -> bowtie2;
      out_2 [label="alignments"];
      bowtie2 -> out_2;
   }    

**Options:**
  - **index** (str, required) -- Path to bowtie2 index (not containing file suffixes).
    
**Required tools:** bowtie2, dd, mkfifo, pigz

**CPU Cores:** 6

.. index:: bowtie2_generate_index

bowtie2_generate_index
======================

bowtie2-build builds a Bowtie index from a set of DNA sequences.
bowtie2-build outputs a set of 6 files with suffixes .1.bt2, .2.bt2, .3.bt2,
.4.bt2, .rev.1.bt2, and .rev.2.bt2. In the case of a large index these
suffixes will have a bt2l termination. These files together constitute the
index: they are all that is needed to align reads to that reference.
The original sequence FASTA files are no longer used by Bowtie 2 once the
index is built.

http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer

typical command line:

.. code-block:: bash

    bowtie2-build [options]* <reference_in> <bt2_index_base>

**Connections:**
  - Input Connection:
    
    - 'in/reference_sequence'
  - Output Connection:
    
    - 'out/bowtie_index'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      bowtie2_generate_index [style=filled, fillcolor="#fce94f"];
      in_0 [label="reference_sequence"];
      in_0 -> bowtie2_generate_index;
      out_1 [label="bowtie_index"];
      bowtie2_generate_index -> out_1;
   }    

**Options:**
  - **bmax** (int, optional) -- The maximum number of suffixes allowed in a block. Allowing more suffixes per block makes indexing faster, but increases peak memory usage. Setting this option overrides any previous setting for --bmax, or --bmaxdivn. Default (in terms of the --bmaxdivn parameter) is --bmaxdivn 4. This is configured automatically by default; use -a/--noauto to configure manually.
    
  - **bmaxdivn** (int, optional) -- The maximum number of suffixes allowed in a block, expressed as a fraction of the length of the reference. Setting this option overrides any previous setting for --bmax, or --bmaxdivn. Default: --bmaxdivn 4. This is configured automatically by default; use -a/--noauto to configure manually.
    
  - **cutoff** (int, optional) -- Index only the first <int> bases of the reference sequences (cumulative across sequences) and ignore the rest.
    
  - **dcv** (int, optional) -- Use <int> as the period for the difference-cover sample. A larger period yields less memory overhead, but may make suffix sorting slower, especially if repeats are present. Must be a power of 2 no greater than 4096. Default: 1024. This is configured automatically by default; use -a/--noauto to configure manually.
    
  - **ftabchars** (int, optional) -- The ftab is the lookup table used to calculate an initial Burrows-Wheeler range with respect to the first <int> characters of the query. A larger <int> yields a larger lookup table but faster query times. The ftab has size 4^(<int>+1) bytes. The default setting is 10 (ftab is 4MB).
    
  - **index-basename** (str, required) -- Base name used for the bowtie2 index.
    
  - **large-index** (bool, optional) -- Force bowtie2-build to build a large index, even if the reference is less than ~ 4 billion nucleotides long.
    
  - **noauto** (bool, optional) -- Disable the default behavior whereby bowtie2-build automatically selects values for the --bmax, --dcv and --packed parameters according to available memory. Instead, user may specify values for those parameters. If memory is exhausted during indexing, an error message will be printed; it is up to the user to try new parameters.
    
  - **nodc** (bool, optional) -- Disable use of the difference-cover sample. Suffix sorting becomes quadratic-time in the worst case (where the worst case is an extremely repetitive reference). Default: off.
    
  - **offrate** (int, optional) -- To map alignments back to positions on the reference sequences, it's necessary to annotate ('mark') some or all of the Burrows-Wheeler rows with their corresponding location on the genome. -o/--offrate governs how many rows get marked: the indexer will mark every 2^<int> rows. Marking more rows makes reference-position lookups faster, but requires more memory to hold the annotations at runtime. The default is 5 (every 32nd row is marked; for human genome, annotations occupy about 340 megabytes).
    
  - **packed** (bool, optional) -- Use a packed (2-bits-per-nucleotide) representation for DNA strings. This saves memory but makes indexing 2-3 times slower. Default: off. This is configured automatically by default; use -a/--noauto to configure manually.
    
  - **seed** (int, optional) -- Use <int> as the seed for pseudo-random number generator.
    
**Required tools:** bowtie2-build, dd, pigz

**CPU Cores:** 6

.. index:: bwa_backtrack

bwa_backtrack
=============

bwa-backtrack is the bwa algorithm designed for Illumina sequence reads up
to 100bp. The computation of the alignments is done by running 'bwa aln'
first, to align the reads, followed by running 'bwa samse' or 'bwa sampe'
afterwards to generate the final SAM output.

http://bio-bwa.sourceforge.net/

typical command line for single-end data:

.. code-block:: bash

    bwa aln <bwa-index> <first-read.fastq> > <first-read.sai>
    bwa samse <bwa-index> <first-read.sai> <first-read.fastq> > <sam-output>

typical command line for paired-end data:

.. code-block:: bash

    bwa aln <bwa-index> <first-read.fastq> > <first-read.sai>
    bwa aln <bwa-index> <second-read.fastq> > <second-read.sai>
    bwa sampe <bwa-index> <first-read.sai> <second-read.sai> <first-read.fastq> <second-read.fastq> > <sam-output>

**Connections:**
  - Input Connection:
    
    - 'in/first_read'
    - 'in/second_read'
  - Output Connection:
    
    - 'out/alignments'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      bwa_backtrack [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> bwa_backtrack;
      in_1 [label="second_read"];
      in_1 -> bwa_backtrack;
      out_2 [label="alignments"];
      bwa_backtrack -> out_2;
   }    

**Options:**
  - **aln-0** (bool, optional) -- When aln-b is specified, only use single-end reads in mapping.
    
  - **aln-1** (bool, optional) -- When aln-b is specified, only use the first read in a read pair in mapping (skip single-end reads and the second reads).
    
  - **aln-2** (bool, optional) -- When aln-b is specified, only use the second read in a read pair in mapping.
    
  - **aln-B** (int, optional) -- Length of barcode starting from the 5'-end. When INT is positive, the barcode of each read will be trimmed before mapping and will be written at the BC SAM tag. For paired-end reads, the barcode from both ends are concatenated. [0]
    
  - **aln-E** (int, optional) -- Gap extension penalty [4]
    
  - **aln-I** (bool, optional) -- The input is in the Illumina 1.3+ read format (quality equals ASCII-64).
    
  - **aln-M** (int, optional) -- Mismatch penalty. BWA will not search for suboptimal hits with a score lower than (bestScore-misMsc). [3]
    
  - **aln-N** (bool, optional) -- Disable iterative search. All hits with no more than maxDiff differences will be found. This mode is much slower than the default.
    
  - **aln-O** (int, optional) -- Gap open penalty [11]
    
  - **aln-R** (int, optional) -- Proceed with suboptimal alignments if there are no more than INT equally best hits. This option only affects paired-end mapping. Increasing this threshold helps to improve the pairing accuracy at the cost of speed, especially for short reads (~32bp).
    
  - **aln-b** (bool, optional) -- Specify the input read sequence file is the BAM format. For paired-end data, two ends in a pair must be grouped together and options aln-1 or aln-2 are usually applied to specify which end should be mapped. Typical command lines for mapping pair-end data in the BAM format are:

    ::

        bwa aln ref.fa -b1 reads.bam > 1.sai
        bwa aln ref.fa -b2 reads.bam > 2.sai
        bwa sampe ref.fa 1.sai 2.sai reads.bam reads.bam > aln.sam

  - **aln-c** (bool, optional) -- Reverse query but not complement it, which is required for alignment in the color space. (Disabled since 0.6.x)
    
  - **aln-d** (int, optional) -- Disallow a long deletion within INT bp towards the 3'-end [16]
    
  - **aln-e** (int, optional) -- Maximum number of gap extensions, -1 for k-difference mode (disallowing long gaps) [-1]
    
  - **aln-i** (int, optional) -- Disallow an indel within INT bp towards the ends [5]
    
  - **aln-k** (int, optional) -- Maximum edit distance in the seed [2]
    
  - **aln-l** (int, optional) -- Take the first INT subsequence as seed. If INT is larger than the query sequence, seeding will be disabled. For long reads, this option is typically ranged from 25 to 35 for '-k 2'. [inf]
    
  - **aln-n** (float, optional) -- Maximum edit distance if the value is INT, or the fraction of missing alignments given 2% uniform base error rate if FLOAT. In the latter case, the maximum edit distance is automatically chosen for different read lengths. [0.04]
    
  - **aln-o** (int, optional) -- Maximum number of gap opens [1]
    
  - **aln-q** (int, optional) -- Parameter for read trimming. BWA trims a read down to argmax_x{\sum_{i=x+1}^l(INT-q_i)} if q_l<INT where l is the original read length. [0]
    
  - **aln-t** (int, optional) -- Number of threads (multi-threading mode) [1]
    
  - **index** (str, required) -- Path to BWA index
    
  - **sampe-N** (int, optional) -- Maximum number of alignments to output in the XA tag for reads paired properly. If a read has more than INT hits, the XA tag will not be written. [3]
    
  - **sampe-P** (bool, optional) -- Load the entire FM-index into memory to reduce disk operations (base-space reads only). With this option, at least 1.25N bytes of memory are required, where N is the length of the genome.
    
  - **sampe-a** (int, optional) -- Maximum insert size for a read pair to be considered being mapped properly. Since 0.4.5, this option is only used when there are not enough good alignment to infer the distribution of insert sizes. [500]
    
  - **sampe-n** (int, optional) -- Maximum number of alignments to output in the XA tag for reads paired properly. If a read has more than INT hits, the XA tag will not be written. [3]
    
  - **sampe-o** (int, optional) -- Maximum occurrences of a read for pairing. A read with more occurrneces will be treated as a single-end read. Reducing this parameter helps faster pairing. [100000]
    
  - **sampe-r** (str, optional) -- Specify the read group in a format like '@RG	ID:foo	SM:bar'. [null]
    
  - **samse-n** (int, optional) -- Maximum number of alignments to output in the XA tag for reads paired properly. If a read has more than INT hits, the XA tag will not be written. [3]
    
  - **samse-r** (str, optional) -- Specify the read group in a format like '@RG	ID:foo	SM:bar'. [null]
    
**Required tools:** bwa, dd, mkfifo, pigz

**CPU Cores:** 8

.. index:: bwa_generate_index

bwa_generate_index
==================

This step generates the index database from sequences in the FASTA format.

Typical command line:

.. code-block:: bash

    bwa index -p <index-basename> <seqeunce.fasta>

**Connections:**
  - Input Connection:
    
    - 'in/reference_sequence'
  - Output Connection:
    
    - 'out/bwa_index'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      bwa_generate_index [style=filled, fillcolor="#fce94f"];
      in_0 [label="reference_sequence"];
      in_0 -> bwa_generate_index;
      out_1 [label="bwa_index"];
      bwa_generate_index -> out_1;
   }    

**Options:**
  - **index-basename** (str, required) -- Prefix of the created index database
    
**Required tools:** bwa, dd, pigz

**CPU Cores:** 6

.. index:: bwa_mem

bwa_mem
=======

Align 70bp-1Mbp query sequences with the BWA-MEM algorithm. Briefly, the
algorithm works by seeding alignments with maximal exact matches (MEMs) and
then extending seeds with the affine-gap Smith-Waterman algorithm (SW).

http://bio-bwa.sourceforge.net/bwa.shtml

Typical command line:

.. code-block:: bash

    bwa mem [options] <bwa-index> <first-read.fastq> [<second-read.fastq>] > <sam-output>

**Connections:**
  - Input Connection:
    
    - 'in/first_read'
    - 'in/second_read'
  - Output Connection:
    
    - 'out/alignments'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      bwa_mem [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> bwa_mem;
      in_1 [label="second_read"];
      in_1 -> bwa_mem;
      out_2 [label="alignments"];
      bwa_mem -> out_2;
   }    

**Options:**
  - **A** (int, optional) -- score for a sequence match, which scales options -TdBOELU unless overridden [1]
    
  - **B** (int, optional) -- penalty for a mismatch [4]
    
  - **C** (bool, optional) -- append FASTA/FASTQ comment to SAM output
    
  - **D** (float, optional) -- drop chains shorter than FLOAT fraction of the longest overlapping chain [0.50]
    
  - **E** (str, optional) -- gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [1,1]
    
  - **H** (str, optional) -- insert STR to header if it starts with @; or insert lines in FILE [null]
    
  - **L** (str, optional) -- penalty for 5'- and 3'-end clipping [5,5]
    
  - **M** (str, optional) -- mark shorter split hits as secondary
    
  - **O** (str, optional) -- gap open penalties for deletions and insertions [6,6]
    
  - **P** (bool, optional) -- skip pairing; mate rescue performed unless -S also in use
    
  - **R** (str, optional) -- read group header line such as '@RG	ID:foo	SM:bar' [null]
    
  - **S** (bool, optional) -- skip mate rescue
    
  - **T** (int, optional) -- minimum score to output [30]
    
  - **U** (int, optional) -- penalty for an unpaired read pair [17]
    
  - **V** (bool, optional) -- output the reference FASTA header in the XR tag
    
  - **W** (int, optional) -- discard a chain if seeded bases shorter than INT [0]
    
  - **Y** (str, optional) -- use soft clipping for supplementary alignments
    
  - **a** (bool, optional) -- output all alignments for SE or unpaired PE
    
  - **c** (int, optional) -- skip seeds with more than INT occurrences [500]
    
  - **d** (int, optional) -- off-diagonal X-dropoff [100]
    
  - **e** (bool, optional) -- discard full-length exact matches
    
  - **h** (str, optional) -- if there are <INT hits with score >80% of the max score, output all in XA [5,200]
    
  - **index** (str, required) -- Path to BWA index
    
  - **j** (bool, optional) -- treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)
    
  - **k** (int, optional) -- minimum seed length [19]
    
  - **m** (int, optional) -- perform at most INT rounds of mate rescues for each read [50]
    
  - **p** (bool, optional) -- smart pairing (ignoring in2.fq)
    
  - **r** (float, optional) -- look for internal seeds inside a seed longer than {-k} * FLOAT [1.5]
    
  - **t** (int, optional) -- number of threads [6]
    
  - **v** (int, optional) -- verbose level: 1=error, 2=warning, 3=message, 4+=debugging [3]
    
  - **w** (int, optional) -- band width for banded alignment [100]
    
  - **x** (str, optional) -- read type. Setting -x changes multiple parameters unless overriden [null]
        pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref)
        ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref)
        intractg: -B9 -O16 -L5  (intra-species contigs to ref)
    
  - **y** (int, optional) -- seed occurrence for the 3rd round seeding [20]
    
**Required tools:** bwa, dd, mkfifo, pigz

**CPU Cores:** 6

.. index:: chromhmm_binarizebam

chromhmm_binarizebam
====================

This command converts coordinates of aligned reads into binarized data form
from which a chromatin state model can be learned. The binarization is based
on a poisson background model. If no control data is specified the parameter
to the poisson distribution is the global average number of reads per bin.
If control data is specified the global average number of reads is
multiplied by the local enrichment for control reads as determined by the
specified parameters. Optionally intermediate signal files can also be
outputted and these signal files can later be directly converted into binary
form using the BinarizeSignal command.

**Connections:**
  - Input Connection:
    
    - 'in/alignments'
  - Output Connection:
    
    - 'out/alignments'
    - 'out/metrics'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      chromhmm_binarizebam [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> chromhmm_binarizebam;
      out_1 [label="alignments"];
      chromhmm_binarizebam -> out_1;
      out_2 [label="metrics"];
      chromhmm_binarizebam -> out_2;
   }    

**Options:**
  - **b** (int, optional)
    
  - **c** (str, optional)
    
  - **center** (bool, optional)
    
  - **chrom_sizes_file** (str, required)
    
  - **control** (dict, required)
    
  - **e** (int, optional)
    
  - **f** (int, optional)
    
  - **g** (int, optional)
    
  - **n** (int, optional)
    
  - **o** (str, optional)
    
  - **p** (float, optional)
    
  - **peaks** (bool, optional)
    
  - **s** (int, optional)
    
  - **strictthresh** (bool, optional)
    
  - **t** (str, optional)
    
  - **u** (int, optional)
    
  - **w** (int, optional)
    
**Required tools:** ChromHMM, echo, ln

**CPU Cores:** 8

.. index:: cutadapt

cutadapt
========

Cutadapt finds and removes adapter sequences, primers, poly-A tails and
other types of unwanted sequence from your high-throughput sequencing reads.

https://cutadapt.readthedocs.org/en/stable/

**Connections:**
  - Input Connection:
    
    - 'in/first_read'
    - 'in/second_read'
  - Output Connection:
    
    - 'out/first_read'
    - 'out/log_first_read'
    - 'out/log_second_read'
    - 'out/second_read'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      cutadapt [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> cutadapt;
      in_1 [label="second_read"];
      in_1 -> cutadapt;
      out_2 [label="first_read"];
      cutadapt -> out_2;
      out_3 [label="log_first_read"];
      cutadapt -> out_3;
      out_4 [label="log_second_read"];
      cutadapt -> out_4;
      out_5 [label="second_read"];
      cutadapt -> out_5;
   }    

**Options:**
  - **adapter-R1** (str, optional) -- Adapter sequence to be clipped off of thefirst read.
    
  - **adapter-R2** (str, optional) -- Adapter sequence to be clipped off of thesecond read
    
  - **adapter-file** (str, optional) -- File containing adapter sequences to be clipped off of the reads.
    
  - **adapter-type** (str, optional)
    
    - possible values: '-a', '-g', '-b'
    
  - **fix_qnames** (bool, required) -- If set to true, only the leftmost string without spaces of the QNAME field of the FASTQ data is kept. This might be necessary for downstream analysis.
    
  - **use_reverse_complement** (bool, required) -- The reverse complement of adapter sequences 'adapter-R1' and 'adapter-R2' are used for adapter clipping.
    
**Required tools:** cat, cutadapt, dd, fix_qnames, mkfifo, pigz

**CPU Cores:** 4

.. index:: fastqc

fastqc
======

The fastqc step  is a wrapper for the fastqc tool. It generates some quality
metrics for fastq files. For this specific instance only the zip archive is
preserved.

http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

**Connections:**
  - Input Connection:
    
    - 'in/first_read'
    - 'in/second_read'
  - Output Connection:
    
    - 'out/first_read_fastqc_report'
    - 'out/first_read_fastqc_report_webpage'
    - 'out/first_read_log_stderr'
    - 'out/second_read_fastqc_report'
    - 'out/second_read_fastqc_report_webpage'
    - 'out/second_read_log_stderr'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      fastqc [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> fastqc;
      in_1 [label="second_read"];
      in_1 -> fastqc;
      out_2 [label="first_read_fastqc_report"];
      fastqc -> out_2;
      out_3 [label="first_read_fastqc_report_webpage"];
      fastqc -> out_3;
      out_4 [label="first_read_log_stderr"];
      fastqc -> out_4;
      out_5 [label="second_read_fastqc_report"];
      fastqc -> out_5;
      out_6 [label="second_read_fastqc_report_webpage"];
      fastqc -> out_6;
      out_7 [label="second_read_log_stderr"];
      fastqc -> out_7;
   }    

**Options:**
**Required tools:** fastqc, mkdir, mv

**CPU Cores:** 1

.. index:: fastx_quality_stats

fastx_quality_stats
===================

fastx_quality_stats generates a text file containing quality information
of the input FASTQ data.

Documentation: |fastx_toolkit_link|

**Connections:**
  - Input Connection:
    
    - 'in/first_read'
    - 'in/second_read'
  - Output Connection:
    
    - 'out/first_read_quality_stats'
    - 'out/second_read_quality_stats'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      fastx_quality_stats [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> fastx_quality_stats;
      in_1 [label="second_read"];
      in_1 -> fastx_quality_stats;
      out_2 [label="first_read_quality_stats"];
      fastx_quality_stats -> out_2;
      out_3 [label="second_read_quality_stats"];
      fastx_quality_stats -> out_3;
   }    

**Options:**
  - **new_output_format** (bool, optional)
    
  - **quality** (int, optional)
    
**Required tools:** cat, dd, fastx_quality_stats, mkfifo, pigz

**CPU Cores:** 4

.. index:: fix_cutadapt

fix_cutadapt
============

This step takes FASTQ data and removes both reads of a paired-end read, if
one of them has been completely removed by cutadapt (or any other software).

**Connections:**
  - Input Connection:
    
    - 'in/first_read'
    - 'in/second_read'
  - Output Connection:
    
    - 'out/first_read'
    - 'out/second_read'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      fix_cutadapt [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> fix_cutadapt;
      in_1 [label="second_read"];
      in_1 -> fix_cutadapt;
      out_2 [label="first_read"];
      fix_cutadapt -> out_2;
      out_3 [label="second_read"];
      fix_cutadapt -> out_3;
   }    

**Options:**
**Required tools:** cat, dd, fix_cutadapt, mkfifo, pigz

**CPU Cores:** 4

.. index:: htseq_count

htseq_count
===========

The htseq-count script counts the number of reads overlapping a feature.
Input needs to be a file with aligned sequencing reads and a list of genomic
features. For more information see: |htseq_link|

**Connections:**
  - Input Connection:
    
    - 'in/alignments'
    - 'in/features'
  - Output Connection:
    
    - 'out/counts'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      htseq_count [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> htseq_count;
      in_1 [label="features"];
      in_1 -> htseq_count;
      out_2 [label="counts"];
      htseq_count -> out_2;
   }    

**Options:**
  - **a** (int, optional)
    
  - **feature-file** (str, optional)
    
  - **idattr** (str, optional)
    
  - **mode** (str, optional)
    
    - possible values: 'union', 'intersection-strict', 'intersection-nonempty'
    
  - **order** (str, required)
    
    - possible values: 'name', 'pos'
    
  - **stranded** (str, required)
    
    - possible values: 'yes', 'no', 'reverse'
    
  - **type** (str, optional)
    
**Required tools:** dd, htseq-count, pigz, samtools

**CPU Cores:** 2

.. index:: macs2

macs2
=====

Model-based Analysis of ChIP-Seq (MACS) is a algorithm, for the identifcation
of transcript factor binding sites. MACS captures the influence of genome
complexity to evaluate the significance of enriched ChIP regions, and MACS
improves the spatial resolution of binding sites through combining the
information of both sequencing tag position and orientation. MACS can be
easily used for ChIP-Seq data alone, or with control sample data to increase
the specificity.

https://github.com/taoliu/MACS

typical command line for single-end data:

.. code-block:: bash

    macs2 callpeak --treatment <aligned-reads> [--control <aligned-reads>] --name <run-id> --gsize 2.7e9

**Connections:**
  - Input Connection:
    
    - 'in/alignments'
  - Output Connection:
    
    - 'out/broadpeaks'
    - 'out/broadpeaks-xls'
    - 'out/diagnosis'
    - 'out/gappedpeaks'
    - 'out/log'
    - 'out/model'
    - 'out/narrowpeaks'
    - 'out/narrowpeaks-xls'
    - 'out/summits'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      macs2 [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> macs2;
      out_1 [label="broadpeaks"];
      macs2 -> out_1;
      out_2 [label="broadpeaks-xls"];
      macs2 -> out_2;
      out_3 [label="diagnosis"];
      macs2 -> out_3;
      out_4 [label="gappedpeaks"];
      macs2 -> out_4;
      out_5 [label="log"];
      macs2 -> out_5;
      out_6 [label="model"];
      macs2 -> out_6;
      out_7 [label="narrowpeaks"];
      macs2 -> out_7;
      out_8 [label="narrowpeaks-xls"];
      macs2 -> out_8;
      out_9 [label="summits"];
      macs2 -> out_9;
   }    

**Options:**
  - **broad** (bool, optional)
    
  - **broad-cutoff** (float, optional)
    
  - **buffer-size** (int, optional)
    
  - **call-summits** (bool, optional)
    
  - **control** (dict, required)
    
  - **down-sample** (bool, optional)
    
  - **format** (str, required)
    
    - possible values: 'AUTO', 'ELAND', 'ELANDMULTI', 'ELANDMULTIPET', 'ELANDEXPORT', 'BED', 'SAM', 'BAM', 'BAMPE', 'BOWTIE'
    
  - **gsize** (str, required)
    
  - **keep-dup** (int, optional)
    
  - **llocal** (str, optional)
    
  - **pvalue** (float, optional)
    
  - **qvalue** (float, optional)
    
  - **read-length** (int, optional)
    
  - **shift** (int, optional)
    
  - **slocal** (str, optional)
    
  - **to-large** (bool, optional)
    
  - **verbose** (int, optional)
    
    - possible values: '0', '1', '2', '3'
    
**Required tools:** macs2, mkdir, mv, pigz

**CPU Cores:** 4

.. index:: merge_fasta_files

merge_fasta_files
=================

This step merges all .fasta(.gz) files belonging to a certain sample.
The output files are gzipped.

**Connections:**
  - Input Connection:
    
    - 'in/sequence'
  - Output Connection:
    
    - 'out/sequence'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      merge_fasta_files [style=filled, fillcolor="#fce94f"];
      in_0 [label="sequence"];
      in_0 -> merge_fasta_files;
      out_1 [label="sequence"];
      merge_fasta_files -> out_1;
   }    

**Options:**
  - **compress-output** (bool, optional) -- If set to true output is gzipped.
    
  - **output-fasta-basename** (str, optional) -- Name used as prefix for FASTA output.
    
**Required tools:** cat, dd, mkfifo, pigz

**CPU Cores:** 12

.. index:: merge_fastq_files

merge_fastq_files
=================

This step merges all .fastq(.gz) files belonging to a certain sample.
First and second read files are merged separately. The output files are
gzipped.

**Connections:**
  - Input Connection:
    
    - 'in/first_read'
    - 'in/second_read'
  - Output Connection:
    
    - 'out/first_read'
    - 'out/second_read'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      merge_fastq_files [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> merge_fastq_files;
      in_1 [label="second_read"];
      in_1 -> merge_fastq_files;
      out_2 [label="first_read"];
      merge_fastq_files -> out_2;
      out_3 [label="second_read"];
      merge_fastq_files -> out_3;
   }    

**Options:**
**Required tools:** cat, dd, mkfifo, pigz

**CPU Cores:** 12

.. index:: narrowpeak_to_bed

narrowpeak_to_bed
=================

**Connections:**
  - Input Connection:
    
    - 'in/broadpeaks'
    - 'in/narrowpeaks'
  - Output Connection:
    
    - 'out/broadpeaks-bed'
    - 'out/narrowpeaks-bed'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      narrowpeak_to_bed [style=filled, fillcolor="#fce94f"];
      in_0 [label="broadpeaks"];
      in_0 -> narrowpeak_to_bed;
      in_1 [label="narrowpeaks"];
      in_1 -> narrowpeak_to_bed;
      out_2 [label="broadpeaks-bed"];
      narrowpeak_to_bed -> out_2;
      out_3 [label="narrowpeaks-bed"];
      narrowpeak_to_bed -> out_3;
   }    

**Options:**
  - **genome** (str, required)
    
  - **sort-by-name** (bool, required)
    
  - **temp-sort-directory** (str, required) -- Intermediate sort files are stored intothis directory.
    
**Required tools:** bedClip, bedtools

**CPU Cores:** 8

.. index:: picard_add_replace_read_groups

picard_add_replace_read_groups
==============================

Replace read groups in a BAM file. This tool enables the user to replace all
read groups in the INPUT file with a single new read group and assign all
reads to this read group in the OUTPUT BAM file.

Documentation: |picard_add_replace_read_groups_link|

**Connections:**
  - Input Connection:
    
    - 'in/alignments'
  - Output Connection:
    
    - 'out/alignments'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      picard_add_replace_read_groups [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> picard_add_replace_read_groups;
      out_1 [label="alignments"];
      picard_add_replace_read_groups -> out_1;
   }    

**Options:**
  - **COMPRESSION_LEVEL** (int, optional) -- Compression level for all compressed files created (e.g. BAM and GELI). Default value: 5. This option can be set to "null" to clear the default value.
    
  - **CREATE_INDEX** (bool, optional) -- Whether to create a BAM index when writing a coordinate-sorted BAM file. Default value: false. This option can be set to "null" to clear the default value. 
    
  - **CREATE_MD5_FILE** (bool, optional) -- Whether to create an MD5 digest for any BAM or FASTQ files created. Default value: false. This option can be set to "null" to clear the default value.
    
  - **GA4GH_CLIENT_SECRETS** (str, optional) -- Google Genomics API client_secrets.json file path. Default value: client_secrets.json. This option can be set to "null" to clear the default value.
    
  - **MAX_RECORDS_IN_RAM** (int, optional) -- When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed. Default value: 500000. This option can be set to "null" to clear the default value.
    
  - **QUIET** (bool, optional) -- Whether to suppress job-summary info on System.err. Default value: false. This option can be set to "null" to clear the default value.
    
  - **REFERENCE_SEQUENCE** (str, optional) -- Reference sequence file. Default value: null.
    
  - **RGCN** (str, optional) -- Read Group sequencing center name. Default value: null.
    
  - **RGDS** (str, optional) -- Read Group description. Default value: null.
    
  - **RGDT** (str, optional) -- Read Group run date. Default value: null.
    
  - **RGID** (str, optional) -- Read Group ID Default value: 1. This option can be set to 'null' to clear the default value.
    
  - **RGLB** (str, required) -- Read Group library
    
  - **RGPG** (str, optional) -- Read Group program group. Default value: null.
    
  - **RGPI** (int, optional) -- Read Group predicted insert size. Default value: null.
    
  - **RGPL** (str, required) -- Read Group platform (e.g. illumina, solid)
    
  - **RGPM** (str, optional) -- Read Group platform model. Default value: null.
    
  - **RGPU** (str, required) -- Read Group platform unit (eg. run barcode)
    
  - **SORT_ORDER** (str, optional) -- Optional sort order to output in. If not supplied OUTPUT is in the same order as INPUT. Default value: null. Possible values: {unsorted, queryname, coordinate, duplicate}
    
    - possible values: 'unsorted', 'queryname', 'coordinate', 'duplicate'
    
  - **TMP_DIR** (str, optional) -- A file. Default value: null. This option may be specified 0 or more times.
    
  - **VALIDATION_STRINGENCY** (str, optional) -- Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded. Default value: STRICT. This option can be set to "null" to clear the default value.
    
    - possible values: 'STRICT', 'LENIENT', 'SILENT'
    
  - **VERBOSITY** (str, optional) -- Control verbosity of logging. Default value: INFO. This option can be set to "null" to clear the default value.
    
    - possible values: 'ERROR', 'WARNING', 'INFO', 'DEBUG'
    
**Required tools:** picard-tools

**CPU Cores:** 6

.. index:: picard_markduplicates

picard_markduplicates
=====================

Identifies duplicate reads.
This tool locates and tags duplicate reads (both PCR and optical/
sequencing-driven) in a BAM or SAM file, where duplicate reads are defined
as originating from the same original fragment of DNA.
Duplicates are identified as read pairs having identical 5' positions
(coordinate and strand) for both reads in a mate pair (and optinally,
matching unique molecular identifier reads; see BARCODE_TAG option).
Optical, or more broadly Sequencing, duplicates are duplicates that appear
clustered together spatially during sequencing and can arise from optical/
imagine-processing artifacts or from bio-chemical processes during clonal
amplification and sequencing; they are identified using the READ_NAME_REGEX
and the OPTICAL_DUPLICATE_PIXEL_DISTANCE options.
The tool's main output is a new SAM or BAM file in which duplicates have
been identified in the SAM flags field, or optionally removed (see
REMOVE_DUPLICATE and REMOVE_SEQUENCING_DUPLICATES), and optionally marked
with a duplicate type in the 'DT' optional attribute.
In addition, it also outputs a metrics file containing the numbers of
READ_PAIRS_EXAMINED, UNMAPPED_READS, UNPAIRED_READS,
UNPAIRED_READ_DUPLICATES, READ_PAIR_DUPLICATES, and
READ_PAIR_OPTICAL_DUPLICATES.

Usage example:

.. code-block:: bash

    java -jar picard.jar MarkDuplicates I=input.bam O=marked_duplicates.bam M=marked_dup_metrics.txt

Documentation: |picard_mark_duplicates_link|

**Connections:**
  - Input Connection:
    
    - 'in/alignments'
  - Output Connection:
    
    - 'out/alignments'
    - 'out/metrics'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      picard_markduplicates [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> picard_markduplicates;
      out_1 [label="alignments"];
      picard_markduplicates -> out_1;
      out_2 [label="metrics"];
      picard_markduplicates -> out_2;
   }    

**Options:**
  - **ASSUME_SORTED** (bool, optional)
    
  - **COMMENT** (str, optional)
    
  - **COMPRESSION_LEVEL** (int, optional) -- Compression level for all compressed files created (e.g. BAM and GELI). Default value: 5. This option can be set to "null" to clear the default value.
    
  - **CREATE_INDEX** (bool, optional) -- Whether to create a BAM index when writing a coordinate-sorted BAM file. Default value: false. This option can be set to "null" to clear the default value. 
    
  - **CREATE_MD5_FILE** (bool, optional) -- Whether to create an MD5 digest for any BAM or FASTQ files created. Default value: false. This option can be set to "null" to clear the default value.
    
  - **GA4GH_CLIENT_SECRETS** (str, optional) -- Google Genomics API client_secrets.json file path. Default value: client_secrets.json. This option can be set to "null" to clear the default value.
    
  - **MAX_FILE_HANDLES** (int, optional)
    
  - **MAX_RECORDS_IN_RAM** (int, optional) -- When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed. Default value: 500000. This option can be set to "null" to clear the default value.
    
  - **OPTICAL_DUPLICATE_PIXEL_DISTANCE** (int, optional)
    
  - **PROGRAM_GROUP_COMMAND_LINE** (str, optional)
    
  - **PROGRAM_GROUP_NAME** (str, optional)
    
  - **PROGRAM_GROUP_VERSION** (str, optional)
    
  - **PROGRAM_RECORD_ID** (str, optional)
    
  - **QUIET** (bool, optional) -- Whether to suppress job-summary info on System.err. Default value: false. This option can be set to "null" to clear the default value.
    
  - **READ_NAME_REGEX** (str, optional)
    
  - **REFERENCE_SEQUENCE** (str, optional) -- Reference sequence file. Default value: null.
    
  - **SORTING_COLLECTION_SIZE_RATIO** (float, optional)
    
  - **TMP_DIR** (str, optional) -- A file. Default value: null. This option may be specified 0 or more times.
    
  - **VALIDATION_STRINGENCY** (str, optional) -- Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded. Default value: STRICT. This option can be set to "null" to clear the default value.
    
    - possible values: 'STRICT', 'LENIENT', 'SILENT'
    
  - **VERBOSITY** (str, optional) -- Control verbosity of logging. Default value: INFO. This option can be set to "null" to clear the default value.
    
    - possible values: 'ERROR', 'WARNING', 'INFO', 'DEBUG'
    
**Required tools:** picard-tools

**CPU Cores:** 12

.. index:: picard_merge_sam_bam_files

picard_merge_sam_bam_files
==========================

Documentation: |picard_merge_sam_files_link|

**Connections:**
  - Input Connection:
    
    - 'in/alignments'
  - Output Connection:
    
    - 'out/alignments'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      picard_merge_sam_bam_files [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> picard_merge_sam_bam_files;
      out_1 [label="alignments"];
      picard_merge_sam_bam_files -> out_1;
   }    

**Options:**
  - **ASSUME_SORTED** (bool, optional) -- If true, assume that the input files are in the same sort order as the requested output sort order, even if their headers say otherwise. Default value: false. This option can be set to 'null' to clear the default value. Possible values: {true, false}
    
  - **COMMENT** (str, optional) -- Comment(s) to include in the merged output file's header. Default value: null.
    
  - **COMPRESSION_LEVEL** (int, optional) -- Compression level for all compressed files created (e.g. BAM and GELI). Default value: 5. This option can be set to "null" to clear the default value.
    
  - **CREATE_INDEX** (bool, optional) -- Whether to create a BAM index when writing a coordinate-sorted BAM file. Default value: false. This option can be set to "null" to clear the default value. 
    
  - **CREATE_MD5_FILE** (bool, optional) -- Whether to create an MD5 digest for any BAM or FASTQ files created. Default value: false. This option can be set to "null" to clear the default value.
    
  - **GA4GH_CLIENT_SECRETS** (str, optional) -- Google Genomics API client_secrets.json file path. Default value: client_secrets.json. This option can be set to "null" to clear the default value.
    
  - **INTERVALS** (str, optional) -- An interval list file that contains the locations of the positions to merge. Assume bam are sorted and indexed. The resulting file will contain alignments that may overlap with genomic regions outside the requested region. Unmapped reads are discarded. Default value: null.
    
  - **MAX_RECORDS_IN_RAM** (int, optional) -- When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed. Default value: 500000. This option can be set to "null" to clear the default value.
    
  - **MERGE_SEQUENCE_DICTIONARIES** (bool, optional) -- Merge the sequence dictionaries. Default value: false. This option can be set to 'null' to clear the default value. Possible values: {true, false}
    
  - **QUIET** (bool, optional) -- Whether to suppress job-summary info on System.err. Default value: false. This option can be set to "null" to clear the default value.
    
  - **REFERENCE_SEQUENCE** (str, optional) -- Reference sequence file. Default value: null.
    
  - **SORT_ORDER** (str, optional) -- Sort order of output file. Default value: coordinate. This option can be set to 'null' to clear the default value. Possible values: {unsorted, queryname, coordinate, duplicate}
    
    - possible values: 'unsorted', 'queryname', 'coordinate', 'duplicate'
    
  - **TMP_DIR** (str, optional) -- A file. Default value: null. This option may be specified 0 or more times.
    
  - **USE_THREADING** (bool, optional) -- Option to create a background thread to encode, compress and write to disk the output file. The threaded version uses about 20% more CPU and decreases runtime by ~20% when writing out a compressed BAM file. Default value: false. This option can be set to 'null' to clear the default value. Possible values: {true, false}
    
  - **VALIDATION_STRINGENCY** (str, optional) -- Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded. Default value: STRICT. This option can be set to "null" to clear the default value.
    
    - possible values: 'STRICT', 'LENIENT', 'SILENT'
    
  - **VERBOSITY** (str, optional) -- Control verbosity of logging. Default value: INFO. This option can be set to "null" to clear the default value.
    
    - possible values: 'ERROR', 'WARNING', 'INFO', 'DEBUG'
    
**Required tools:** ln, picard-tools

**CPU Cores:** 12

.. index:: preseq_complexity_curve

preseq_complexity_curve
=======================

The preseq package is aimed at predicting the yield of distinct reads from a
genomic library from an initial sequencing experiment. The estimates can then
be used to examine the utility of further sequencing, optimize the sequencing
depth, or to screen multiple libraries to avoid low complexity samples.

c_curve computes the expected yield of distinct reads for experiments smaller
than the input experiment in a .bed or .bam file through resampling. The full
set of parameters can be outputed by simply typing the program name. If
output.txt is the desired output file name and input.bed is the input .bed
file, then simply type:

.. code-block:: bash

    preseq c_curve -o output.txt input.sort.bed

**Connections:**
  - Input Connection:
    
    - 'in/alignments'
  - Output Connection:
    
    - 'out/complexity_curve'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      preseq_complexity_curve [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> preseq_complexity_curve;
      out_1 [label="complexity_curve"];
      preseq_complexity_curve -> out_1;
   }    

**Options:**
  - **hist** (bool, optional) -- input is a text file containing the observed histogram
    
  - **pe** (bool, required) -- input is paired end read file
    
  - **seg_len** (int, optional) -- maximum segment length when merging paired end bam reads (default: 5000)
    
  - **step** (int, optional) -- step size gin extrapolations (default: 1e+06)
    
  - **vals** (bool, optional) -- input is a text file containing only the observed counts
    
**Required tools:** preseq

**CPU Cores:** 4

.. index:: preseq_future_genome_coverage

preseq_future_genome_coverage
=============================

The preseq package is aimed at predicting the yield of distinct reads from a
genomic library from an initial sequencing experiment. The estimates can then
be used to examine the utility of further sequencing, optimize the sequencing
depth, or to screen multiple libraries to avoid low complexity samples.

gc_extrap computes the expected genomic coverage for deeper sequencing for
single cell sequencing experiments. The input should be a mr or bed file.
The tool bam2mr is provided to convert sorted bam or sam files to mapped
read format.

**Connections:**
  - Input Connection:
    
    - 'in/alignments'
  - Output Connection:
    
    - 'out/future_genome_coverage'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      preseq_future_genome_coverage [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> preseq_future_genome_coverage;
      out_1 [label="future_genome_coverage"];
      preseq_future_genome_coverage -> out_1;
   }    

**Options:**
  - **bin_size** (int, optional) -- bin size (default: 10)
    
  - **bootstraps** (int, optional) -- number of bootstraps (default: 100)
    
  - **cval** (float, optional) -- level for confidence intervals (default: 0.95)
    
  - **extrap** (int, optional) -- maximum extrapolation in base pairs (default: 1e+12)
    
  - **max_width** (int, optional) -- max fragment length, set equal to read length for single end reads
    
  - **quick** (bool, optional) -- quick mode: run gc_extrap without bootstrapping for confidence intervals
    
  - **step** (int, optional) -- step size in bases between extrapolations (default: 1e+08)
    
  - **terms** (int, optional) -- maximum number of terms
    
**Required tools:** preseq

**CPU Cores:** 4

.. index:: preseq_future_yield

preseq_future_yield
===================

The preseq package is aimed at predicting the yield of distinct reads from a
genomic library from an initial sequencing experiment. The estimates can then
be used to examine the utility of further sequencing, optimize the sequencing
depth, or to screen multiple libraries to avoid low complexity samples.

lc_extrap computes the expected future yield of distinct reads and bounds on
the number of total distinct reads in the library and the associated
confidence intervals.

**Connections:**
  - Input Connection:
    
    - 'in/alignments'
  - Output Connection:
    
    - 'out/future_yield'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      preseq_future_yield [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> preseq_future_yield;
      out_1 [label="future_yield"];
      preseq_future_yield -> out_1;
   }    

**Options:**
  - **bootstraps** (int, optional) -- number of bootstraps (default: 100)
    
  - **cval** (float, optional) -- level for confidence intervals (default: 0.95)
    
  - **dupl_level** (float, optional) -- fraction of duplicate to predict (default: 0.5)
    
  - **extrap** (int, optional) -- maximum extrapolation (default: 1e+10)
    
  - **hist** (bool, optional) -- input is a text file containing the observed histogram
    
  - **pe** (bool, required) -- input is paired end read file
    
  - **quick** (bool, optional) -- quick mode, estimate yield without bootstrapping for confidence intervals
    
  - **seg_len** (int, optional) -- maximum segment length when merging paired end bam reads (default: 5000)
    
  - **step** (int, optional) -- step size in extrapolations (default: 1e+06)
    
  - **terms** (int, optional) -- maximum number of terms
    
  - **vals** (bool, optional) -- input is a text file containing only the observed counts
    
**Required tools:** preseq

**CPU Cores:** 4

.. index:: remove_duplicate_reads_runs

remove_duplicate_reads_runs
===========================

Duplicates are removed by Picard tools 'MarkDuplicates'.

typical command line:

.. code-block:: bash

    MarkDuplicates INPUT=<SAM/BAM> OUTPUT=<SAM/BAM> METRICS_FILE=<metrics-out> REMOVE_DUPLICATES=true

**Connections:**
  - Input Connection:
    
    - 'in/alignments'
  - Output Connection:
    
    - 'out/alignments'
    - 'out/metrics'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      remove_duplicate_reads_runs [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> remove_duplicate_reads_runs;
      out_1 [label="alignments"];
      remove_duplicate_reads_runs -> out_1;
      out_2 [label="metrics"];
      remove_duplicate_reads_runs -> out_2;
   }    

**Options:**
**Required tools:** MarkDuplicates

**CPU Cores:** 12

.. index:: rseqc

rseqc
=====

The RSeQC step can be used to evaluate aligned reads in a BAM file. RSeQC
does not only report raw sequence-based metrics, but also quality control
metrics like read distribution, gene coverage, and sequencing depth.

**Connections:**
  - Input Connection:
    
    - 'in/alignments'
  - Output Connection:
    
    - 'out/bam_stat'
    - 'out/infer_experiment'
    - 'out/read_distribution'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      rseqc [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> rseqc;
      out_1 [label="bam_stat"];
      rseqc -> out_1;
      out_2 [label="infer_experiment"];
      rseqc -> out_2;
      out_3 [label="read_distribution"];
      rseqc -> out_3;
   }    

**Options:**
  - **reference** (str, required) -- Reference gene model in bed fomat. [required]
    
**Required tools:** bam_stat.py, cat, infer_experiment.py, read_distribution.py

**CPU Cores:** 1

.. index:: sam_to_sorted_bam

sam_to_sorted_bam
=================

The step sam_to_sorted_bam builds on 'samtools sort' to sort SAM files and
output BAM files.

Sort alignments by leftmost coordinates, or by read name when -n is used.
An appropriate @HD-SO sort order header tag will be added or an existing
one updated if necessary.

Documentation: |samtools_link|

**Connections:**
  - Input Connection:
    
    - 'in/alignments'
  - Output Connection:
    
    - 'out/alignments'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      sam_to_sorted_bam [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> sam_to_sorted_bam;
      out_1 [label="alignments"];
      sam_to_sorted_bam -> out_1;
   }    

**Options:**
  - **genome-faidx** (str, required)
    
  - **sort-by-name** (bool, required)
    
  - **temp-sort-directory** (str, required) -- Intermediate sort files are stored intothis directory.
    
**Required tools:** dd, pigz, samtools

**CPU Cores:** 8

.. index:: samtools_faidx

samtools_faidx
==============

Index reference sequence in the FASTA format or extract subsequence from
indexed reference sequence. If no region is specified, faidx will index the
file and create <ref.fasta>.fai on the disk. If regions are specified, the
subsequences will be retrieved and printed to stdout in the FASTA format.

**Connections:**
  - Input Connection:
    
    - 'in/sequence'
  - Output Connection:
    
    - 'out/indices'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      samtools_faidx [style=filled, fillcolor="#fce94f"];
      in_0 [label="sequence"];
      in_0 -> samtools_faidx;
      out_1 [label="indices"];
      samtools_faidx -> out_1;
   }    

**Options:**
**Required tools:** mv, samtools

**CPU Cores:** 4

.. index:: samtools_index

samtools_index
==============

Index a coordinate-sorted BAM or CRAM file for fast random access.
(Note that this does not work with SAM files even if they are bgzip
compressed to index such files, use tabix(1) instead.)

Documentation: |samtools_link|

**Connections:**
  - Input Connection:
    
    - 'in/alignments'
  - Output Connection:
    
    - 'out/alignments'
    - 'out/index_stats'
    - 'out/indices'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      samtools_index [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> samtools_index;
      out_1 [label="alignments"];
      samtools_index -> out_1;
      out_2 [label="index_stats"];
      samtools_index -> out_2;
      out_3 [label="indices"];
      samtools_index -> out_3;
   }    

**Options:**
  - **index_type** (str, required)
    
    - possible values: 'bai', 'csi'
    
**Required tools:** ln, samtools

**CPU Cores:** 4

.. index:: samtools_stats

samtools_stats
==============

samtools stats collects statistics from BAM files and outputs in a text
format. The output can be visualized graphically using plot-bamstats.

Documentation: |samtools_link|

**Connections:**
  - Input Connection:
    
    - 'in/alignments'
  - Output Connection:
    
    - 'out/stats'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      samtools_stats [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> samtools_stats;
      out_1 [label="stats"];
      samtools_stats -> out_1;
   }    

**Options:**
**Required tools:** dd, pigz, samtools

**CPU Cores:** 1

.. index:: segemehl

segemehl
========

segemehl is a software to map short sequencer reads to reference genomes.
Unlike other methods, segemehl is able to detect not only mismatches but
also insertions and deletions. Furthermore, segemehl is not limited to a
specific read length and is able to mapprimer- or polyadenylation
contaminated reads correctly.

This step creates at first two FIFOs. The first is used to provide the
genome data for segemehl and the second is used for the output of the
unmapped reads:

.. code-block:: bash

    mkfifo genome_fifo unmapped_fifo
    cat <genome-fasta> -o genome_fifo

The executed segemehl command is this:

.. code-block:: bash

    segemehl -d genome_fifo -i <genome-index-file> -q <read1-fastq> [-p <read2-fastq>] -u unmapped_fifo -H 1 -t 11 -s -S -D 0 -o /dev/stdout |  pigz --blocksize 4096 --processes 2 -c

The unmapped reads are saved via these commands:

.. code-block:: bash

    cat unmapped_fifo | pigz --blocksize 4096 --processes 2 -c > <unmapped-fastq>

**Connections:**
  - Input Connection:
    
    - 'in/first_read'
    - 'in/second_read'
  - Output Connection:
    
    - 'out/alignments'
    - 'out/log'
    - 'out/unmapped'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      segemehl [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> segemehl;
      in_1 [label="second_read"];
      in_1 -> segemehl;
      out_2 [label="alignments"];
      segemehl -> out_2;
      out_3 [label="log"];
      segemehl -> out_3;
      out_4 [label="unmapped"];
      segemehl -> out_4;
   }    

**Options:**
  - **MEOP** (bool, optional) -- output MEOP field for easier variance calling in SAM (XE:Z:)
    
  - **SEGEMEHL** (bool, optional) -- output SEGEMEHL format (needs to be selected for brief)
    
  - **accuracy** (int, optional) -- min percentage of matches per read in semi-global alignment (default:90)
    
  - **autoclip** (bool, optional) -- autoclip unknown 3prime adapter
    
  - **bisulfite** (int, optional) -- bisulfite mapping with methylC-seq/Lister et al. (=1) or bs-seq/Cokus et al. protocol (=2) (default:0)
    
    - possible values: '0', '1', '2'
    
  - **brief** (bool, optional) -- brief output
    
  - **clipacc** (int, optional) -- clipping accuracy (default:70)
    
  - **differences** (int, optional) -- search seeds initially with <n> differences (default:1)
    
  - **dropoff** (int, optional) -- dropoff parameter for extension (default:8)
    
  - **evalue** (float, optional) -- max evalue (default:5.000000)
    
  - **extensionpenalty** (int, optional) -- penalty for a mismatch during extension (default:4)
    
  - **extensionscore** (int, optional) -- score of a match during extension (default:2)
    
  - **fix-qnames** (bool, optional) -- The QNAMES field of the input will be purged from spaces and everything thereafter.
    
  - **genome** (str, required) -- Path to genome file
    
  - **hardclip** (bool, optional) -- enable hard clipping
    
  - **hitstrategy** (int, optional) -- report only best scoring hits (=1) or all (=0) (default:1)
    
    - possible values: '0', '1'
    
  - **index** (str, required) -- Path to genome index for segemehl
    
  - **jump** (int, optional) -- search seeds with jump size <n> (0=automatic) (default:0)
    
  - **maxinsertsize** (int, optional) -- maximum size of the inserts (paired end) (default:5000)
    
  - **maxinterval** (int, optional) -- maximum width of a suffix array interval, i.e. a query seed will be omitted if it matches more than <n> times (default:100)
    
  - **maxsplitevalue** (float, optional) -- max evalue for splits (default:50.000000)
    
  - **minfraglen** (int, optional) -- min length of a spliced fragment (default:20)
    
  - **minfragscore** (int, optional) -- min score of a spliced fragment (default:18)
    
  - **minsize** (int, optional) -- minimum size of queries (default:12)
    
  - **minsplicecover** (int, optional) -- min coverage for spliced transcripts (default:80)
    
  - **nohead** (bool, optional) -- do not output header
    
  - **order** (bool, optional) -- sorts the output by chromsome and position (might take a while!)
    
  - **polyA** (bool, optional) -- clip polyA tail
    
  - **prime3** (str, optional) -- add 3' adapter (default:none)
    
  - **prime5** (str, optional) -- add 5' adapter (default:none)
    
  - **showalign** (bool, optional) -- show alignments
    
  - **silent** (bool, optional) -- shut up!
    
  - **splicescorescale** (float, optional) -- report spliced alignment with score s only if <f>*s is larger than next best spliced alignment (default:1.000000)
    
  - **splits** (bool, optional) -- detect split/spliced reads (default:none)
    
**Required tools:** cat, dd, fix_qnames, mkfifo, pigz, segemehl

**CPU Cores:** 12

.. index:: segemehl_generate_index

segemehl_generate_index
=======================

The step segemehl_generate_index generates a index for given reference
sequences.

Documentation: |segemehl_link|

**Connections:**
  - Input Connection:
    
    - 'in/reference_sequence'
  - Output Connection:
    
    - 'out/log'
    - 'out/segemehl_index'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      segemehl_generate_index [style=filled, fillcolor="#fce94f"];
      in_0 [label="reference_sequence"];
      in_0 -> segemehl_generate_index;
      out_1 [label="log"];
      segemehl_generate_index -> out_1;
      out_2 [label="segemehl_index"];
      segemehl_generate_index -> out_2;
   }    

**Options:**
  - **index-basename** (str, required) -- Basename for created segemehl index.
    
**Required tools:** dd, mkfifo, pigz, segemehl

**CPU Cores:** 4

.. index:: tophat2

tophat2
=======

TopHat is a fast splice junction mapper for RNA-Seq reads.
It aligns RNA-Seq reads to mammalian-sized genomes using the ultra
high-throughput short read aligner Bowtie, and then analyzes the mapping
results to identify splice junctions between exons.

http://tophat.cbcb.umd.edu/

typical command line:

.. code-block:: bash

    tophat [options]* <index_base> <reads1_1[,...,readsN_1]> [reads1_2,...readsN_2]

**Connections:**
  - Input Connection:
    
    - 'in/first_read'
    - 'in/second_read'
  - Output Connection:
    
    - 'out/align_summary'
    - 'out/alignments'
    - 'out/deletions'
    - 'out/insertions'
    - 'out/junctions'
    - 'out/log_stderr'
    - 'out/misc_logs'
    - 'out/prep_reads'
    - 'out/unmapped'

.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      tophat2 [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> tophat2;
      in_1 [label="second_read"];
      in_1 -> tophat2;
      out_2 [label="align_summary"];
      tophat2 -> out_2;
      out_3 [label="alignments"];
      tophat2 -> out_3;
      out_4 [label="deletions"];
      tophat2 -> out_4;
      out_5 [label="insertions"];
      tophat2 -> out_5;
      out_6 [label="junctions"];
      tophat2 -> out_6;
      out_7 [label="log_stderr"];
      tophat2 -> out_7;
      out_8 [label="misc_logs"];
      tophat2 -> out_8;
      out_9 [label="prep_reads"];
      tophat2 -> out_9;
      out_10 [label="unmapped"];
      tophat2 -> out_10;
   }    

**Options:**
  - **index** (str, required) -- Path to genome index for tophat2
    
  - **library_type** (str, required) -- The default is unstranded (fr-unstranded). If either fr-firststrand or fr-secondstrand is specified, every read alignment will have an XS attribute tag as explained below. Consider supplying library type options below to select the correct RNA-seq protocol.(https://ccb.jhu.edu/software/tophat/manual.shtml)
    
    - possible values: 'fr-unstranded', 'fr-firststrand', 'fr-secondstrand'
    
**Required tools:** mkdir, mv, tar, tophat2

**CPU Cores:** 6

.. |bowtie2_link| raw:: html

   <a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml" target="_blank">http://bowtie-bio.sourceforge.net/bowtie2/index.shtml</a>

.. |fastx_toolkit_link| raw:: html

   <a href="http://hannonlab.cshl.edu/fastx_toolkit/" target="_blank">http://hannonlab.cshl.edu/fastx_toolkit/</a>

.. |htseq_link| raw:: html

   <a href="http://www-huber.embl.de/users/anders/HTSeq/doc/count.html" target="_blank">http://www-huber.embl.de/users/anders/HTSeq/doc/count.html</a>

.. |picard_add_replace_read_groups_link| raw:: html

   <a href="https://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups" target="_blank">https://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups</a>

.. |picard_mark_duplicates_link| raw:: html

   <a href="https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates" target="_blank">https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates</a>

.. |picard_merge_sam_files_link| raw:: html

   <a href="https://broadinstitute.github.io/picard/command-line-overview.html#MergeSamFiles" target="_blank">https://broadinstitute.github.io/picard/command-line-overview.html#MergeSamFiles</a>

.. |samtools_link| raw:: html

   <a href="http://www.htslib.org/doc/samtools.html" target="_blank">http://www.htslib.org/doc/samtools.html</a>

.. |segemehl_link| raw:: html

   <a href="http://www.bioinf.uni-leipzig.de/Software/segemehl/" target="_blank">http://www.bioinf.uni-leipzig.de/Software/segemehl</a>
