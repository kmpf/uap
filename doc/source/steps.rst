###############
Available steps
###############

************
Source steps
************

.. index:: bcl2fastq2_source

bcl2fastq2_source
=================


**Output Connection**
  - **out/bcl2fastq2_log_stderr**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      bcl2fastq2_source [style=filled, fillcolor="#fce94f"];
      out_0 [label="bcl2fastq2_log_stderr"];
      bcl2fastq2_source -> out_0;
   }

**Options:**
  - **adapter-stringency** (float, optional) -- adapter stringency (=0.9)

  - **aggregated-tiles** (str, optional) -- tiles aggregation flag determining structure of input files (=AUTO). recognized values: AUTO - Try to detect correct setting. YES - Tiles are aggregated into single input file. NO - There are separate input files for individual tiles

  - **barcode-mismatches** (str, optional) -- number of allowed mismatches per index multiple entries, comma delimited entries, allowed; each entry is applied to the corresponding index;last entry applies to all remaining indices

  - **create-fastq-for-index-reads** (bool, optional) -- create FASTQ files also for index reads

  - **fastq-compression-level** (int, optional) -- Zlib compression level (1-9) used for FASTQ files (=4)

  - **find-adapters-with-sliding-window** (bool, optional) -- find adapters with simple sliding window algorithm

  - **ignore-missing-bcls** (bool, optional) -- assume N for missing calls

  - **ignore-missing-controls** (bool, optional) -- assume 0 for missing controls

  - **ignore-missing-filter** (bool, optional) -- assume true for missing filters

  - **ignore-missing-positions** (bool, optional) -- assume [0,i] for missing positions, where i is incremented starting from 0

  - **input-dir** (str, optional) -- path to input directory (=<runfolder-dir>/Data/Intensities/BaseCalls/)

  - **intensities-dir** (str, optional) -- path to intensities directory (=<input-dir>/../). If intensities directory is specified, also input directory must be specified.

  - **interop-dir** (str, optional) -- path to demultiplexing statistics directory (=<runfolder-dir>/InterOp/)

  - **loading-threads** (int, optional) -- number of threads used for loading BCL data (=4)

  - **mask-short-adapter-reads** (int, optional) -- smallest number of remaining bases (after masking bases below the minimum trimmed read length) below which whole read is masked (=22)

  - **min-log-level** (str, optional) -- minimum log level recognized values: NONE, FATAL, ERROR, WARNING, INFO, DEBUG, TRACE. (INFO)

  - **minimum-trimmed-read-length** (int, optional) -- minimum read length after adapter trimming (=35)

  - **no-bgzf-compression** (bool, optional) -- Turn off BGZF compression for FASTQ files

  - **no-lane-splitting** (bool, optional) -- Do not split fastq files by lane.

  - **output-dir** (str, required)
  - **processing-threads** (int, optional) -- number of threads used for processing demultipled data (=100% of available CPUs)

  - **reports-dir** (str, optional) -- path to reporting directory (=<output-dir>/Reports/)

  - **runfolder-dir** (str, required) -- path to runfolder directory (=./)

  - **sample-sheet** (str, optional) -- path to the sample sheet(=<runfolder-dir>/SampleSheet.csv)

  - **stats-dir** (str, optional) -- path to human-readable demultiplexing statistics directory (=<runfolder-dir>/InterOp/)

  - **tiles** (str, optional) -- Comma-separated list of regular expressions to select only a subset of the tiles available in the flow-cell.Multiple entries allowed, each applies to the corresponding base-calls.For example: * to select all the tiles ending with 5 in all lanes: tiles [0-9][0-9][0-9]5. * to select tile 2 in lane 1 and all the tiles in the other lanes: tiles s_1_0002,s_[2-8]

  - **use-bases-mask** (str, optional) -- Specifies how to use each cycle.

  - **with-failed-reads** (bool, optional) -- include non-PF clusters

  - **write-fastq-reverse-complement** (bool, optional) -- Generate FASTQs containing reverse complements of actual data

  - **writing-threads** (int, optional) -- number of threads used for writing FASTQ data (=4)


**Required tools:** bcl2fastq, mkdir (coreutils)

This step provides input files which already exists and therefore creates no tasks in the pipeline.

.. index:: bcl2fastq_source

bcl2fastq_source
================


**Output Connection**
  - **out/make_log_stderr**
  - **out/configureBcl2Fastq_log_stderr**
  - **out/sample_sheet**


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

**Required tools:** configureBclToFastq.pl, make, mkdir (coreutils), mv (coreutils)

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

**Output Connection**
  - **out/second_read** (optional)
  - **out/first_read**


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
      out_1 [label="second_read", style=filled, fillcolor="#a7a7a7"];
      fastq_source -> out_1;
   }

**Options:**
  - **first_read** (str, required) -- Part of the file name that marks all files containing sequencing data of the first read. Example: 'R1.fastq' or '_1.fastq'

  - **group** (str, optional) -- A regular expression which is applied to found files, and which is used to determine the sample name from the file name. For example, ``(Sample_\d+)_R[12].fastq.gz``, when applied to a file called ``Sample_1_R1.fastq.gz``, would result in a sample name of ``Sample_1``. You can specify multiple capture groups in the regular expression.

  - **indices** (str/dict, optional) -- path to a CSV file or a dictionary of sample_id: barcode entries.

  - **paired_end** (bool, optional) -- Specify whether the samples are paired end or not.

  - **pattern** (str, optional) -- A file name pattern, for example ``/home/test/fastq/Sample_*.fastq.gz``.

  - **sample_id_prefix** (str, optional) -- This optional prefix is prepended to every sample name.

  - **sample_to_files_map** (dict/str, optional) -- A listing of sample names and their associated files. This must be provided as a YAML dictionary.

  - **second_read** (str, optional) -- Part of the file name that marks all files containing sequencing data of the second read. Example: 'R2.fastq' or '_2.fastq'


This step provides input files which already exists and therefore creates no tasks in the pipeline.

.. index:: fetch_chrom_sizes_source

fetch_chrom_sizes_source
========================


**Output Connection**
  - **out/chromosome_sizes**


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


**Required tools:** cp (coreutils), fetchChromSizes

This step provides input files which already exists and therefore creates no tasks in the pipeline.

.. index:: raw_file_source

raw_file_source
===============


**Output Connection**
  - **out/raw**


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

  - **sample_id_prefix** (str, optional) -- This optional prefix is prepended to every sample name.

  - **sample_to_files_map** (dict/str, optional) -- A listing of sample names and their associated files. This must be provided as a YAML dictionary.


This step provides input files which already exists and therefore creates no tasks in the pipeline.

.. index:: raw_file_sources

raw_file_sources
================


**Output Connection**
  - **out/raws**


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
  - **group** (str, required) -- **This is a LEGACY step.** Do NOT use it, better use the ``raw_file_source`` step. A regular expression which is applied to found files, and which is used to determine the sample name from the file name. For example, ``(Sample_\d+)_R[12].fastq.gz``, when applied to a file called ``Sample_1_R1.fastq.gz``, would result in a sample name of ``Sample_1``. You can specify multiple capture groups in the regular expression.

  - **paired_end** (bool, required) -- Specify whether the samples are paired end or not.

  - **pattern** (str, required) -- A file name pattern, for example ``/home/test/fastq/Sample_*.fastq.gz``.

  - **sample_id_prefix** (str, optional) -- This optional prefix is prepended to every sample name.


This step provides input files which already exists and therefore creates no tasks in the pipeline.

.. index:: raw_url_source

raw_url_source
==============


**Output Connection**
  - **out/raw**


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
  - **dd-blocksize** (str, optional)    - default value: 256k

  - **filename** (str, optional) -- local file name of downloaded file

  - **hashing-algorithm** (str, optional) -- hashing algorithm to use
    - possible values: 'md5', 'sha1', 'sha224', 'sha256', 'sha384', 'sha512'


  - **path** (str, required) -- directory to move downloaded file to

  - **secure-hash** (str, optional) -- expected secure hash of downloaded file

  - **uncompress** (bool, optional) -- File is uncompressed after download

  - **url** (str, required) -- Download URL


**Required tools:** cp (coreutils), curl, dd (coreutils), mkdir (coreutils), pigz

This step provides input files which already exists and therefore creates no tasks in the pipeline.

.. index:: raw_url_sources

raw_url_sources
===============


**Output Connection**
  - **out/raw**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      raw_url_sources [style=filled, fillcolor="#fce94f"];
      out_0 [label="raw"];
      raw_url_sources -> out_0;
   }

**Options:**
  - **dd-blocksize** (str, optional)    - default value: 256k

  - **run-download-info** (dict, required) -- Dictionary of dictionaries. The keys are the names of the runs. The values are dictionaries whose keys are identical with the options of an 'raw_url_source' source step. An example: <name>: filename: <filename> hashing-algorithm: <hashing-algorithm> path: <path> secure-hash: <secure-hash> uncompress: <uncompress> url: <url>


**Required tools:** cp (coreutils), curl, dd (coreutils), pigz

This step provides input files which already exists and therefore creates no tasks in the pipeline.

.. index:: run_folder_source

run_folder_source
=================



    This source looks for fastq.gz files in
    ``[path]/Unaligned/Project_*/Sample_*`` and pulls additional information
    from CSV sample sheets it finds. It also makes sure that index information
    for all samples is coherent and unambiguous.

**Output Connection**
  - **out/second_read** (optional)
  - **out/first_read**


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
      out_1 [label="second_read", style=filled, fillcolor="#a7a7a7"];
      run_folder_source -> out_1;
   }

**Options:**
  - **first_read** (str, required) -- Part of the file name that marks all files containing sequencing data of the first read. Example: '_R1.fastq' or '_1.fastq'
    - default value: _R1

  - **paired_end** (bool, optional) -- Is the project a paired-end sequencing project?

  - **path** (str, required) -- Path to the sequencing directories that contain the fastq[.gz] files.

  - **project** (str, optional) -- Name of the project. If provided, this is appendedto the path string
    - default value: *

  - **project_name** (str, optional) -- Name of the project. If provided, this is appendedto the path string. This option has the same meaning as 'project',however, the prefix 'Project\_' is not added. If 'project' and 'project_name' are provided, 'project_name' is choosen.
    - default value: *

  - **samples** (str, optional) -- Pattern for the sample directory names inside path/[Project\_]project[_name]
    - default value: Sample_*

  - **second_read** (str, optional) -- Part of the file name that marks all files containing sequencing data of the second read. Example: 'R2.fastq.gz' or '_2.fastq'
    - default value: _R2

  - **unaligned_included** (bool, optional) -- Is the typical Unaligned folder included in path?
    - default value: True


This step provides input files which already exists and therefore creates no tasks in the pipeline.

****************
Processing steps
****************

.. index:: adapterremoval

adapterremoval
==============



    AdapterRemoval (ver. 2.1.7)

    This program searches for and removes remnant adapter sequences from
    your read data.  The program can analyze both single end and paired end
    data.  For detailed explanation of the parameters, please refer to the
    man page.  For comments, suggestions  and feedback please contact Stinus
    Lindgreen (stinus@binf.ku.dk) and Mikkel Schubert (MikkelSch@gmail.com).

    If you use the program, please cite the paper:
    Schubert, Lindgreen, and Orlando (2016). AdapterRemoval v2: rapid
    adapter trimming, identification, and read merging.
    BMC Research Notes, 12;9(1):88.

    http://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-016-1900-2

    "Pipeline specific "input and output expected to be gzipped"

**Input Connection**
  - **in/second_read** (optional)
  - **in/first_read**

**Output Connection**
  - **out/settings**
  - **out/collapsed** (optional)
  - **out/discarded**
  - **out/log_stdout**
  - **out/log_stderr**
  - **out/truncated** (optional) Format: **fastq** - Truncated single end reads.
  - **out/pair1.truncated** (optional) Format: **fastq** - Truncated first read of paired end reads.
  - **out/singleton.truncated** (optional)
  - **out/pair2.truncated** (optional) Format: **fastq** - Truncated first secind of paired end reads.
  - **out/collapsed.truncated** (optional)


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      adapterremoval [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> adapterremoval;
      in_1 [label="second_read", style=filled, fillcolor="#a7a7a7"];
      in_1 -> adapterremoval;
      out_2 [label="collapsed", style=filled, fillcolor="#a7a7a7"];
      adapterremoval -> out_2;
      out_3 [label="collapsed.truncated", style=filled, fillcolor="#a7a7a7"];
      adapterremoval -> out_3;
      out_4 [label="discarded"];
      adapterremoval -> out_4;
      out_5 [label="log_stderr"];
      adapterremoval -> out_5;
      out_6 [label="log_stdout"];
      adapterremoval -> out_6;
      out_7 [label="pair1.truncated", style=filled, fillcolor="#a7a7a7"];
      adapterremoval -> out_7;
      out_8 [label="pair2.truncated", style=filled, fillcolor="#a7a7a7"];
      adapterremoval -> out_8;
      out_9 [label="settings"];
      adapterremoval -> out_9;
      out_10 [label="singleton.truncated", style=filled, fillcolor="#a7a7a7"];
      adapterremoval -> out_10;
      out_11 [label="truncated", style=filled, fillcolor="#a7a7a7"];
      adapterremoval -> out_11;
   }

**Options:**
  - **adapter1** (str, required) -- Adapter sequence expected to be found in mate 1 reads [current: AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG]

  - **adapter2** (str, optional) -- Adapter sequence expected to be found in mate 2 reads [current: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT]

  - **collapse** (bool, required) -- When set, paired ended read alignments of --minalignmentlength or more bases are combined into a single consensus sequence, representing the complete insert, and written to either basename.collapsed or basename.collapsed.truncated (if trimmed due to low-quality bases following collapse); for single-ended reads, putative complete inserts are identified as having at least --minalignmentlength bases overlap with the adapter sequence, and are written to the the same files [current:off]

  - **cores** (int, optional) -- workaround to specify cores for grid engine and threads ie
    - default value: 1

  - **identify-adapters** (bool, optional)
  - **maxlength** (int, optional) -- Reads longer than this length are discarded following trimming [current: 4294967295]

  - **maxns** (int, optional) -- Reads containing more ambiguous bases (N) than this number after trimming are discarded [current: 1000]

  - **minadapteroverlap** (bool, optional) -- In single-end mode, reads are only trimmed if the overlap between read and the adapter is at least X bases long, not counting ambiguous nucleotides (N); this is independant of the --minalignmentlength when using --collapse, allowing a conservative selection of putative complete inserts while ensuring that all possible adapter contamination is trimmed [current: 0].

  - **minalignmentlength** (bool, optional) -- If --collapse is set, paired reads must overlap at least this number of bases to be collapsed, and single-ended reads must overlap at least this number of bases with the adapter to be considered complete template molecules [current: 11]

  - **minlength** (int, optional) -- Reads shorter than this length are discarded following trimming [current: 15]

  - **minquality** (int, optional) -- Inclusive minimum; see --trimqualities for details [current: 2]

  - **mm** (int, optional) -- Max error-rate when aligning reads and/or adapters. If > 1, the max error-rate is set to 1 / MISMATCH_RATE; if < 0, the defaults are used, otherwise the user-supplied value is used directly. [defaults: 1/3 for trimming; 1/10 when identifing adapters]

  - **qualitybase** (str, optional) -- Quality base used to encode Phred scores in input; either 33, 64, or solexa [current: 33]
    - possible values: '33', '64', 'solexa'


  - **seed** (int, optional)    - default value: 22595

  - **shift** (int, optional) -- Consider alignments where up to N nucleotides are missing from the 5' termini [current: 2]

  - **trimns** (bool, optional) -- If set, trim ambiguous bases (N) at 5'/3' termini [current: off]

  - **trimqualities** (bool, optional) -- If set, trim bases at 5'/3' termini with quality scores <= to --minquality value [current: off]


**Required tools:** adapterremoval, mv (coreutils)

**CPU Cores:** 1

.. index:: bam_to_bedgraph_and_bigwig

bam_to_bedgraph_and_bigwig
==========================


**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/bigwig**
  - **out/bedgraph**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      bam_to_bedgraph_and_bigwig [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> bam_to_bedgraph_and_bigwig;
      out_1 [label="bedgraph"];
      bam_to_bedgraph_and_bigwig -> out_1;
      out_2 [label="bigwig"];
      bam_to_bedgraph_and_bigwig -> out_2;
   }

**Options:**
  - **chromosome-sizes** (str, required)
  - **temp-sort-dir** (str, optional)

**Required tools:** bedGraphToBigWig, bedtools, sort (coreutils)

**CPU Cores:** 8

.. index:: bam_to_genome_browser

bam_to_genome_browser
=====================


**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/alignments**


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
  - **bedtools-genomecov-split** (bool, required)    - default value: True

  - **bedtools-genomecov-strand** (str, optional)    - possible values: '+', '-'


  - **chromosome-sizes** (str, required)
  - **dd-blocksize** (str, optional)    - default value: 256k

  - **output-format** (str, required)    - default value: bigWig
    - possible values: 'bed', 'bigBed', 'bedGraph', 'bigWig'


  - **trackline** (dict, optional)
  - **trackopts** (dict, optional)

**Required tools:** bedGraphToBigWig, bedToBigBed, bedtools, dd (coreutils), mkfifo (coreutils), pigz

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

    The input reads must come as .f[ast]q[.gz] files.

    http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

    typical command line::

        bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r>} -S [<hit>]

    This step wraps release: bowtie2

**Input Connection**
  - **in/second_read** (optional)
  - **in/first_read**

**Output Connection**
  - **out/unaligned** (optional) -  unpaired reads that didn't align
  - **out/un-conc** (optional) - pairs that didn't align concordantly
  - **out/al-conc** (optional) - pairs that aligned concordantly at least once
  - **out/met-file** (optional) - metrics file
  - **out/alignments**
  - **out/log_stderr**
  - **out/al** (optional) - unpaired reads that aligned at least once


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
      in_1 [label="second_read", style=filled, fillcolor="#a7a7a7"];
      in_1 -> bowtie2;
      out_2 [label="al", style=filled, fillcolor="#a7a7a7"];
      bowtie2 -> out_2;
      out_3 [label="al-conc", style=filled, fillcolor="#a7a7a7"];
      bowtie2 -> out_3;
      out_4 [label="alignments"];
      bowtie2 -> out_4;
      out_5 [label="log_stderr"];
      bowtie2 -> out_5;
      out_6 [label="met-file", style=filled, fillcolor="#a7a7a7"];
      bowtie2 -> out_6;
      out_7 [label="un-conc", style=filled, fillcolor="#a7a7a7"];
      bowtie2 -> out_7;
      out_8 [label="unaligned", style=filled, fillcolor="#a7a7a7"];
      bowtie2 -> out_8;
   }

**Options:**
  - **D** (int, optional) -- give up extending after <int> failed extends in a row (default=15)

  - **L** (int, optional) -- length of seed substrings; must be >3, <32 (default=22)

  - **N** (int, optional) -- Max # mismatches in seed alignment; can be 0 or 1 (default=0)

  - **R** (int, optional) -- for reads w/ repetitive seeds, try <int> sets of seeds (default=2)

  - **al** (str, optional) -- Write unpaired reads that aligned at least once to connection out/al

  - **al-conc** (str, optional) -- write pairs that aligned concordantly at least once to out/al-conc

  - **all** (bool, optional) -- report all alignments; very slow, MAPQ not meaningful

  - **compress** (bool, optional) -- Use pigz to compress bowtie2 results.
    - default value: True

  - **cores** (int, optional) -- number of alignment threads to launch (default=1)

  - **dd-blocksize** (str, optional)    - default value: 2M

  - **dpad** (int, optional) -- include <int> extra ref chars on sides of DP table (default=15)

  - **end-to-end** (bool, optional) -- entire read must align; no clipping. Decide forthis or local option. (default)

  - **fast** (bool, optional) -- Preset, same as: -D 10 -R 2 -N 0 -L 22 -i S,0,2.50

  - **fast-local** (bool, optional) -- Preset, same as: -D 10 -R 2 -N 0 -L 22 -i S,1,1.75

  - **ff** (bool, optional) -- -1, -2 mates align fw/fw

  - **fifo** (bool, optional) -- Use dd and pigz to pipe into bowtie2 with fifos. This does not work reliably with bowtie <= 2.3.4.2 due to a race condition (http://seqanswers.com/forums/showthread.php?t=16540).

  - **fr** (bool, optional) -- -1, -2 mates align fw/rev (default)

  - **gbar** (int, optional) -- disallow gaps within <int> nucs of read extremes (default=4)

  - **i** (str, optional) -- interval between seed substrings w/r/t read len (default="S,1,1.15")

  - **ignore-quals** (bool, optional) -- treat all quality values as 30 on Phred scale

  - **index** (str, required) -- Path to bowtie2 index (not containing file suffixes).

  - **int-quals** (bool, optional) -- Qualities encoded as space-delimited integers

  - **k** (int, optional) -- report up to <int> alns per read; MAPQ not meaningful

  - **local** (bool, optional) -- local alignment; ends might be soft clipped. Decide for this or end-to-end option.

  - **ma** (int, optional) -- match bonus (0 for end-to-end, 2 for local). (default=0)

  - **maxins** (int, optional) -- maximum fragment length (default=500)

  - **met** (int, optional) -- report internal counters & metrics every <int> secs (default=1)

  - **met-file** (bool, optional) -- send metrics to file to connection out/met-file

  - **met-stderr** (bool, optional) -- send metrics to stderr

  - **minins** (int, optional) -- minimum fragment length (default=0)

  - **mm** (bool, optional) -- use memory-mapped I/O for index; many 'bowtie's can share

  - **mp** (int, optional) -- max penalty for mismatch; lower qual = lower penalty (default=6)

  - **n-ceil** (str, optional) -- func for max # non-A/C/G/Ts permitted in aln (default="L,0,0.15")

  - **no-1mm-upfront** (bool, optional) -- do not allow 1 mismatch alignments before attempting to scan for the optimal seeded alignments

  - **no-contain** (bool, optional) -- not concordant when one mate alignment contains other

  - **no-discordant** (bool, optional) -- suppress discordant alignments for paired reads

  - **no-dovetail** (bool, optional) -- not concordant when mates extend past each other

  - **no-head** (bool, optional) -- suppress header lines, i.e. lines starting with @

  - **no-mixed** (bool, optional) -- suppress unpaired alignments for paired reads

  - **no-overlap** (bool, optional) -- not concordant when mates overlap at all

  - **no-sq** (bool, optional) -- suppress @SQ header lines

  - **no-unal** (bool, optional) -- suppress SAM records for unaligned reads

  - **nofw** (bool, optional) -- do not align forward (original) version of read

  - **non-deterministic** (bool, optional) -- seed rand. gen. arbitrarily instead of using read attributes

  - **norc** (bool, optional) -- do not align reverse-complement version of read

  - **np** (int, optional) -- penalty for non-A/C/G/Ts in read/ref (default=1)

  - **omit-sec-seq** (bool, optional) -- put * in SEQ and QUAL fields for secondary alignments

  - **phred33** (bool, optional) -- Qualities are Phred+33 (default)

  - **phred64** (bool, optional) -- Qualities are Phred+64

  - **pigz-blocksize** (str, optional)    - default value: 2048

  - **qc-filter** (bool, optional) -- filter out reads that are bad according to QSEQ filter

  - **quiet** (bool, optional) -- print nothing to stderr except serious errors

  - **rdg** (str, optional) -- read gap open, extend penalties (default="5,3")

  - **reorder** (bool, optional) -- force SAM output order to match order of input reads

  - **rf** (bool, optional) -- -1, -2 mates align rev/fw

  - **rfg** (str, optional) -- reference gap open, extend penalties(default="5,3")

  - **rg** (str, optional) -- add <text> (lab:value) to @RG line of SAM header.Note: @RG line only printed when --rg-id is set.

  - **rg-id** (str, optional) -- set read group id, reflected in @RG line and RG:Z: opt field

  - **score-min** (str, optional) -- min acceptable alignment score w/r/t read length(G,20,8 for local, L,-0.6,-0.6 for end-to-end)(default="L,-0.6,-0.6")

  - **seed** (int, optional) -- seed for random number generator (default=0)

  - **sensitive** (bool, optional) -- Preset, same as: -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 (default)

  - **sensitive-local** (bool, optional) -- Preset, same as: -D 15 -R 2 -N 0 -L 20 -i S,1,0.75 (default)

  - **skip** (int, optional) -- Skip the first <int> reads/pairs in the input. Default: none

  - **time** (bool, optional) -- print wall-clock time taken by search phases

  - **trim3** (int, optional) -- Trim <int> bases from 3'/right end of reads(default=0)

  - **trim5** (int, optional) -- Trim <int> bases from 5'/left end of reads (default=0)

  - **un-conc** (str, optional) -- Write pairs that didn't align concordantly to connection out/un-conc

  - **unaligned** (bool, optional) -- Write unpaired reads that didn't align to connection out/unaligned

  - **upto** (int, optional) -- Stop after the first <int> reads/pairs in the input. Default: no limit.

  - **very-fast** (bool, optional) -- Preset, same as: -D 5 -R 1 -N 0 -L 22 -i S,0,2.50

  - **very-fast-local** (bool, optional) -- Preset, same as: -D 5 -R 1 -N 0 -L 25 -i S,1,2.00

  - **very-sensitive** (bool, optional) -- Preset, same as: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50

  - **very-sensitive-local** (bool, optional) -- Preset, same as: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50


**Required tools:** bowtie2, dd (coreutils), mkfifo (coreutils), pigz

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

    typical command line::

        bowtie2-build [options]* <reference_in> <bt2_index_base>

**Input Connection**
  - **in/reference_sequence**

**Output Connection**
  - **out/bowtie_index**


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

  - **dd-blocksize** (str, optional)    - default value: 2M

  - **ftabchars** (int, optional) -- The ftab is the lookup table used to calculate an initial Burrows-Wheeler range with respect to the first <int> characters of the query. A larger <int> yields a larger lookup table but faster query times. The ftab has size 4^(<int>+1) bytes. The default setting is 10 (ftab is 4MB).

  - **index-basename** (str, required) -- Base name used for the bowtie2 index.

  - **large-index** (bool, optional) -- Force bowtie2-build to build a large index, even if the reference is less than 4 billion nucleotides long.

  - **noauto** (bool, optional) -- Disable the default behavior whereby bowtie2-build automatically selects values for the --bmax, --dcv and --packed parameters according to available memory. Instead, user may specify values for those parameters. If memory is exhausted during indexing, an error message will be printed; it is up to the user to try new parameters.

  - **nodc** (bool, optional) -- Disable use of the difference-cover sample. Suffix sorting becomes quadratic-time in the worst case (where the worst case is an extremely repetitive reference). Default: off.

  - **offrate** (int, optional) -- To map alignments back to positions on the reference sequences, it's necessary to annotate ('mark') some or all of the Burrows-Wheeler rows with their corresponding location on the genome. -o/--offrate governs how many rows get marked: the indexer will mark every 2^<int> rows. Marking more rows makes reference-position lookups faster, but requires more memory to hold the annotations at runtime. The default is 5 (every 32nd row is marked; for human genome, annotations occupy about 340 megabytes).

  - **packed** (bool, optional) -- Use a packed (2-bits-per-nucleotide) representation for DNA strings. This saves memory but makes indexing 2-3 times slower. Default: off. This is configured automatically by default; use -a/--noauto to configure manually.

  - **pigz-blocksize** (str, optional)    - default value: 2048

  - **seed** (int, optional) -- Use <int> as the seed for pseudo-random number generator.

  - **threads** (int, optional) -- By default bowtie2-build is using only one thread. Increasing the number of threads will speed up the index building considerably in most cases.


**Required tools:** bowtie2-build, dd (coreutils), pigz

**CPU Cores:** 6

.. index:: bwa_backtrack

bwa_backtrack
=============



    bwa-backtrack is the bwa algorithm designed for Illumina sequence reads up
    to 100bp. The computation of the alignments is done by running 'bwa aln'
    first, to align the reads, followed by running 'bwa samse' or 'bwa sampe'
    afterwards to generate the final SAM output.

    http://bio-bwa.sourceforge.net/

    typical command line for single-end data::

        bwa aln <bwa-index> <first-read.fastq> > <first-read.sai>
        bwa samse <bwa-index> <first-read.sai> <first-read.fastq> > <sam-output>

    typical command line for paired-end data::

        bwa aln <bwa-index> <first-read.fastq> > <first-read.sai>
        bwa aln <bwa-index> <second-read.fastq> > <second-read.sai>
        bwa sampe <bwa-index> <first-read.sai> <second-read.sai>                   <first-read.fastq> <second-read.fastq> > <sam-output>


**Input Connection**
  - **in/second_read**
  - **in/first_read**

**Output Connection**
  - **out/alignments**


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

  - **aln-b** (bool, optional) -- Specify the input read sequence file is the BAM format. For paired-end data, two ends in a pair must be grouped together and options aln-1 or aln-2 are usually applied to specify which end should be mapped. Typical command lines for mapping pair-end data in the BAM format are: bwa aln ref.fa -b1 reads.bam > 1.sai bwa aln ref.fa -b2 reads.bam > 2.sai bwa sampe ref.fa 1.sai 2.sai reads.bam reads.bam > aln.sam 

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
    - default value: 1

  - **dd-blocksize** (str, optional)    - default value: 2M

  - **index** (str, required) -- Path to BWA index

  - **pigz-blocksize** (str, optional)    - default value: 2048

  - **sampe-N** (int, optional) -- Maximum number of alignments to output in the XA tag for disconcordant read pairs (excluding singletons). If a read has more than INT hits, the XA tag will not be written. [10]

  - **sampe-P** (bool, optional) -- Load the entire FM-index into memory to reduce disk operations (base-space reads only). With this option, at least 1.25N bytes of memory are required, where N is the length of the genome.

  - **sampe-a** (int, optional) -- Maximum insert size for a read pair to be considered being mapped properly. Since 0.4.5, this option is only used when there are not enough good alignment to infer the distribution of insert sizes. [500]

  - **sampe-n** (int, optional) -- Maximum number of alignments to output in the XA tag for reads paired properly. If a read has more than INT hits, the XA tag will not be written. [3]

  - **sampe-o** (int, optional) -- Maximum occurrences of a read for pairing. A read with more occurrneces will be treated as a single-end read. Reducing this parameter helps faster pairing. [100000]

  - **sampe-r** (str, optional) -- Specify the read group in a format like '@RG ID:foo SM:bar'. [null]

  - **samse-n** (int, optional) -- Maximum number of alignments to output in the XA tag for reads paired properly. If a read has more than INT hits, the XA tag will not be written. [3]

  - **samse-r** (str, optional) -- Specify the read group in a format like '@RG ID:foo SM:bar'. [null]


**Required tools:** bwa, dd (coreutils), mkfifo (coreutils), pigz

**CPU Cores:** 8

.. index:: bwa_generate_index

bwa_generate_index
==================



    This step generates the index database from sequences in the FASTA format.

    Typical command line::

        bwa index -p <index-basename> <seqeunce.fasta>

**Input Connection**
  - **in/reference_sequence**

**Output Connection**
  - **out/bwa_index**


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


**Required tools:** bwa

**CPU Cores:** 6

.. index:: bwa_mem

bwa_mem
=======



    Align 70bp-1Mbp query sequences with the BWA-MEM algorithm. Briefly, the
    algorithm works by seeding alignments with maximal exact matches (MEMs) and
    then extending seeds with the affine-gap Smith-Waterman algorithm (SW).

    http://bio-bwa.sourceforge.net/bwa.shtml

    Typical command line::

        bwa mem [options] <bwa-index> <first-read.fastq> [<second-read.fastq>]         > <sam-output>

**Input Connection**
  - **in/second_read**
  - **in/first_read**

**Output Connection**
  - **out/alignments**


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

  - **R** (str, optional) -- read group header line such as '@RG ID:foo SM:bar' [null]

  - **S** (bool, optional) -- skip mate rescue

  - **T** (int, optional) -- minimum score to output [30]

  - **U** (int, optional) -- penalty for an unpaired read pair [17]

  - **V** (bool, optional) -- output the reference FASTA header in the XR tag

  - **W** (int, optional) -- discard a chain if seeded bases shorter than INT [0]

  - **Y** (str, optional) -- use soft clipping for supplementary alignments

  - **a** (bool, optional) -- output all alignments for SE or unpaired PE

  - **c** (int, optional) -- skip seeds with more than INT occurrences [500]

  - **d** (int, optional) -- off-diagonal X-dropoff [100]

  - **dd-blocksize** (str, optional)    - default value: 256k

  - **e** (bool, optional) -- discard full-length exact matches

  - **h** (str, optional) -- if there are <INT hits with score >80% of the max score, output all in XA [5,200]

  - **index** (str, required) -- Path to BWA index

  - **j** (bool, optional) -- treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)

  - **k** (int, optional) -- minimum seed length [19]

  - **m** (int, optional) -- perform at most INT rounds of mate rescues for each read [50]

  - **p** (bool, optional) -- smart pairing (ignoring in2.fq)

  - **r** (float, optional) -- look for internal seeds inside a seed longer than {-k} * FLOAT [1.5]

  - **t** (int, optional) -- number of threads [6]
    - default value: 6

  - **v** (int, optional) -- verbose level: 1=error, 2=warning, 3=message, 4+=debugging [3]

  - **w** (int, optional) -- band width for banded alignment [100]

  - **x** (str, optional) -- read type. Setting -x changes multiple parameters unless overriden [null]:: pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0 (PacBio reads to ref) ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0 (Oxford Nanopore 2D-reads to ref) intractg: -B9 -O16 -L5 (intra-species contigs to ref)

  - **y** (int, optional) -- seed occurrence for the 3rd round seeding [20]


**Required tools:** bwa, dd (coreutils), mkfifo (coreutils), pigz

**CPU Cores:** 6

.. index:: cat_text

cat_text
========



    cats text files together

**Input Connection**
  - **in/text**

**Output Connection**
  - **out/text**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      cat_text [style=filled, fillcolor="#fce94f"];
      in_0 [label="text"];
      in_0 -> cat_text;
      out_1 [label="text"];
      cat_text -> out_1;
   }

**Options:**
  - **additionalFiles** (list, optional)
  - **filenameEnding** (str, required)
  - **run_id** (str, optional)    - default value: merged


**Required tools:** cat (coreutils)

**CPU Cores:** 8

.. index:: chimpipe

chimpipe
========




    ChimPipe is a tool to discover gene fusions in human paired-end RNA-Seq data.

    Paper: https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-2-r12

    Paper:
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5209911/

    Manual including typical usage:
    http://chimpipe.readthedocs.io/en/latest/index.html

**Input Connection**
  - **in/second_read**
  - **in/first_read**

**Output Connection**
  - **out/tar_archive**
  - **out/log_stdout**
  - **out/log_stderr**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      chimpipe [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> chimpipe;
      in_1 [label="second_read"];
      in_1 -> chimpipe;
      out_2 [label="log_stderr"];
      chimpipe -> out_2;
      out_3 [label="log_stdout"];
      chimpipe -> out_3;
      out_4 [label="tar_archive"];
      chimpipe -> out_4;
   }

**Options:**
  - **annotation** (str, required) -- Reference gene annotation in .gtf

  - **consensus_seq** (str, optional) -- Sequence pair of consensus splice site bases

  - **cores** (str, required)    - default value: 6

  - **genome_index** (str, required) -- Reference genome index in .gem

  - **library_type** (str, optional) -- Type of sequence library

  - **sample_ID** (str, required) -- Identifier used in output file names

  - **similarity** (str, optional) -- Path to gene pair similarity file

  - **transcriptome_index** (str, required) -- Annotated transcriptome index in .gem

  - **transcriptome_keys** (str, required) -- Transcriptome to genome conversion keys


**Required tools:** chimpipe, mkdir (coreutils), rm (coreutils), tar

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

**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/metrics**
  - **out/alignments**


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
  - **b** (int, optional) -- The number of base pairs in a bin determining the resolution of the model learning and segmentation. By default this parameter value is set to 200 base pairs.

  - **c** (str, optional) -- A directory containing the control input files. If this is not specified then the inputbamdir is used. If no control files are specified then by default a uniform background will be used in determining the binarization thresholds.

  - **center** (bool, optional) -- If this flag is present then the center of the interval is used to determine the bin to assign a read. This can make sense to use if the coordinates are based on already extended reads. If this option is selected, then the strand information of a read and the shift parameter are ignored. By default reads are assigned to a bin based on the position of its 5' end as determined from the strand of the read after shifting an amount determined by the -n shift option.

  - **chrom_sizes_file** (str, required) -- File containing chromosome size information generated by 'fetchChromSizes'

  - **control** (dict, required)
  - **e** (int, optional) -- Specifies the amount that should be subtracted from the end coordinate of a read so that both coordinates are inclusive and 0 based. The default value is 1 corresponding to standard bed convention of the end interval being 0-based but not inclusive.

  - **f** (int, optional) -- This indicates a threshold for the fold enrichment over expected that must be met or exceeded by the observed count in a bin for a present call. The expectation is determined in the same way as the mean parameter for the poission distribution in terms of being based on a uniform background unless control data is specified. This parameter can be useful when dealing with very deeply and/or unevenly sequenced data. By default this parameter value is 0 meaning effectively it is not used.

  - **g** (int, optional) -- This indicates a threshold for the signal that must be met or exceeded by the observed count in a bin for a present call. This parameter can be useful when desiring to directly place a threshold on the signal. By default this parameter value is 0 meaning effectively it is not used.

  - **n** (int, optional) -- The number of bases a read should be shifted to determine a bin assignment. Bin assignment is based on the 5' end of a read shifted this amount with respect to the strand orientation. By default this value is 100.

  - **o** (str, optional) -- This specifies the directory to which control data should be printed. The files will be named CELL_CHROM_controlsignal.txt. Control data will only be outputted if there are control bed files present and an output control directory is specified.

  - **p** (float, optional) -- This option specifies the tail probability of the poisson distribution that the binarization threshold should correspond to. The default value of this parameter is 0.0001.

  - **peaks** (bool, optional) -- This option specifies to treat the bed files as peak calls directly and give a '1' call to any bin overlapping a peak call.

  - **s** (int, optional) -- The amount that should be subtracted from the interval start coordinate so the interval is inclusive and 0 based. Default is 0 corresponding to the standard bed convention.

  - **strictthresh** (bool, optional) -- If this flag is present then the poisson threshold must be strictly greater than the tail probability, otherwise by default the largest integer count for which the tail includes the poisson threshold probability is used.

  - **t** (str, optional)
  - **u** (int, optional) -- An integer pseudocount that is uniformly added to every bin in the control data in order to smooth the control data from 0. The default value is 1.

  - **w** (int, optional) -- This determines the extent of the spatial smoothing in computing the local enrichment for control reads. The local enrichment for control signal in the x-th bin on the chromosome after adding pseudocountcontrol is computed based on the average control counts for all bins within x-w and x+w. If no controldir is specified, then this option is ignored. The default value is 5.


**Required tools:** ChromHMM, echo, ln (coreutils)

**CPU Cores:** 8

.. index:: chromhmm_learnmodel

chromhmm_learnmodel
===================



    This command takes a directory with a set of binarized data files and learns
    a chromatin state model. Binarized data files have "_binary" in the file
    name. The format for the binarized data files are that the first line
    contains the name of the cell separated by a tab with the name of the
    chromosome. The second line contains in tab delimited form the name of each
    mark. The remaining lines correspond to consecutive bins on the chromosome.
    The remaining lines in tab delimited form corresponding to each mark, with a
    "1" for a present call or "0" for an absent call and a "2" if the data is
    considered missing at that interval for the mark.

**Input Connection**
  - **in/chromhmm_binarization**
  - **in/cellmarkfiletable**

**Output Connection**
  - **out/chromhmm_model**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      chromhmm_learnmodel [style=filled, fillcolor="#fce94f"];
      in_0 [label="cellmarkfiletable"];
      in_0 -> chromhmm_learnmodel;
      in_1 [label="chromhmm_binarization"];
      in_1 -> chromhmm_learnmodel;
      out_2 [label="chromhmm_model"];
      chromhmm_learnmodel -> out_2;
   }

**Options:**
  - **assembly** (str, required) -- specifies the genome assembly. overlap and neighborhood enrichments will be called with default parameters using this genome assembly.Assembly names are e.g. hg18, hg19, GRCh38

  - **b** (int, optional) -- The number of base pairs in a bin determining the resolution of the model learning and segmentation. By default this parameter value is set to 200 base pairs.

  - **color** (str, optional) -- This specifies the color of the heat map. "r,g,b" are integer values between 0 and 255 separated by commas. By default this parameter value is 0,0,255 corresponding to blue.

  - **d** (float, optional) -- The threshold on the change on the estimated log likelihood that if it falls below this value, then parameter training will terminate. If this value is less than 0 then it is not used as part of the stopping criteria. The default value for this parameter is 0.001.

  - **e** (float, optional) -- This parameter is only applicable if the load option is selected for the init parameter. This parameter controls the smoothing away from 0 when loading a model. The emission value used in the model initialization is a weighted average of the value in the file and a uniform probability over the two possible emissions. The value in the file gets weight (1-loadsmoothemission) while uniform gets weight loadsmoothemission. The default value of this parameter is 0.02.

  - **h** (float, optional) -- A smoothing constant away from 0 for all parameters in the information based initialization. This option is ignored if random or load are selected for the initialization method. The default value of this parameter is 0.02.

  - **holdcolumnorder** (bool, optional) -- Including this flag suppresses the reordering of the mark columns in the emission parameter table display.

  - **init** (str, optional) -- This specifies the method for parameter initialization method. 'information' is the default method described in (Ernst and Kellis, Nature Methods 2012). 'random' - randomly initializes the parameters from a uniform distribution. 'load' loads the parameters specified in '-m modelinitialfile' and smooths them based on the value of the 'loadsmoothemission' and 'loadsmoothtransition' parameters. The default is information.
    - possible values: 'information', 'random', 'load'


  - **l** (str, optional) -- This file specifies the length of the chromosomes. It is a two column tab delimited file with the first column specifying the chromosome name and the second column the length. If this file is provided then no end coordinate will exceed what is specified in this file. By default BinarizeBed excludes the last partial bin along the chromosome, but if that is included in the binarized data input files then this file should be included to give a valid end coordinate for the last interval.

  - **m** (str, optional) -- This specifies the model file containing the initial parameters which can then be used with the load option

  - **nobed** (bool, optional) -- If this flag is present, then this suppresses the printing of segmentation information in the four column format. The default is to generate a four column segmentation file

  - **nobrowser** (bool, optional) -- If this flag is present, then browser files are not printed. If -nobed is requested then browserfile writing is also suppressed.

  - **noenrich** (bool, optional) -- If this flag is present, then enrichment files are not printed. If -nobed is requested then enrichment file writing is also suppressed.

  - **numstates** (int, required)
  - **r** (int, optional) -- This option specifies the maximum number of iterations over all the input data in the training. By default this is set to 200.

  - **s** (int, optional) -- This allows the specification of the random seed. Randomization is used to determine the visit order of chromosomes in the incremental expectation-maximization algorithm used to train the parameters and also used to generate the initial values of the parameters if random is specified for the init method.

  - **stateordering** (str, optional) -- This determines whether the states are ordered based on the emission or transition parameters. See (Ernst and Kellis, Nature Methods) for details. Default is 'emission'.
    - possible values: 'emission', 'transition'


  - **t** (float, optional) -- This parameter is only applicable if the load option is selected for the init parameter. This parameter controls the smoothing away from 0 when loading a model. The transition value used in the model initialization is a weighted average of the value in the file and a uniform probability over the transitions. The value in the file gets weight (1-loadsmoothtransition) while uniform gets weight loadsmoothtransition. The default value is 0.5.

  - **x** (int, optional) -- This parameter specifies the maximum number of seconds that can be spent optimizing the model parameters. If it is less than 0, then there is no limit and termination is based on maximum number of iterations or a log likelihood change criteria. The default value of this parameter is -1.

  - **z** (int, optional) -- This parameter determines the threshold at which to set extremely low transition probabilities to 0 durining training. Setting extremely low transition probabilities makes model learning more efficient with essentially no impact on the final results. If a transition probability falls below 10^-zerotransitionpower during training it is set to 0. Making this parameter to low and thus the cutoff too high can potentially cause some numerical instability. By default this parameter is set to 8.


**Required tools:** ChromHMM, ls (coreutils), mkdir (coreutils), rm (coreutils), tar, xargs

**CPU Cores:** 8

.. index:: collect_scs

collect_scs
===========



    supa custom succesive aligment info gather thingy

**Input Connection**
  - **in/scs_metrics**

**Output Connection**
  - **out/yaml**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      collect_scs [style=filled, fillcolor="#fce94f"];
      in_0 [label="scs_metrics"];
      in_0 -> collect_scs;
      out_1 [label="yaml"];
      collect_scs -> out_1;
   }

**Options:**
  - **library-type** (str, required)    - possible values: 'unstranded', 'firststranded', 'secondstranded'


  - **rrna-aln-pos** (str, required)
  - **types** (list, required)

**CPU Cores:** 1

.. index:: copy_file

copy_file
=========



    copies a file or a list of files defined by there
    dependencies and filenames

**Input Connection**
  - **in/sequence**

**Output Connection**
  - **out/copied**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      copy_file [style=filled, fillcolor="#fce94f"];
      in_0 [label="sequence"];
      in_0 -> copy_file;
      out_1 [label="copied"];
      copy_file -> out_1;
   }

**Required tools:** cp (coreutils)

**CPU Cores:** 1

.. index:: count_rRNA

count_rRNA
==========



    Tbla bla bla    http://www.htslib.org/doc/samtools.html

**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/report_rRNA**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      count_rRNA [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> count_rRNA;
      out_1 [label="report_rRNA"];
      count_rRNA -> out_1;
   }

**Required tools:** cut (coreutils), grep, pigz, samtools, sort (coreutils), uniq (coreutils)

**CPU Cores:** 8

.. index:: cuffcompare

cuffcompare
===========



    CuffCompare is part of the 'Cufflinks suite of tools' for
    differential expr. analysis of RNA-Seq data and their
    visualisation. This step compares a cufflinks assembly to
    known annotation. Cuffcompare provides classification,
    reference annotation mapping and various statistics for
    Cufflinks transfrags. For details about cuffcompare we refer
    to the author's webpage:

    http://cole-trapnell-lab.github.io/cufflinks/cuffcompare/


**Input Connection**
  - **in/features**

**Output Connection**
  - **out/loci**
  - **out/features**
  - **out/stats**
  - **out/log_stderr**
  - **out/tracking**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      cuffcompare [style=filled, fillcolor="#fce94f"];
      in_0 [label="features"];
      in_0 -> cuffcompare;
      out_1 [label="features"];
      cuffcompare -> out_1;
      out_2 [label="loci"];
      cuffcompare -> out_2;
      out_3 [label="log_stderr"];
      cuffcompare -> out_3;
      out_4 [label="stats"];
      cuffcompare -> out_4;
      out_5 [label="tracking"];
      cuffcompare -> out_5;
   }

**Options:**
  - **C** (bool, optional) -- Enables the "contained" transcripts to be also written in the .combined.gtffile, with the attribute "contained_in" showing the first container transfrag found. By default, without this option, cuffcompare does not write in that file isoforms that were found to be fully contained/covered (with the same compatible intron structure) by other transfrags in the same locus.transcripts

  - **F** (bool, optional) -- Do not discard intron-redundant transfrags if they share the 5p end (if they differ only at the 3p end)

  - **G** (bool, optional) -- generic GFF input file(s): do not assume Cufflinks GTF, do not discard any intron-redundant transfrags

  - **M** (bool, optional) -- Discard (ignore) single-exon transfrags and reference transcripts

  - **N** (bool, optional) -- Discard (ignore) single-exon reference 

  - **Q** (bool, optional) -- For "-r" option, consider only the input transcripts that overlap any of the reference transcripts (Sp-correction)

  - **R** (bool, optional) -- For "-r" option, consider only the reference transcripts that overlap any of the input transfrags (Sn-correction)

  - **V** (bool, optional) -- verbose processing mode (showing all GFF parsing warnings)

  - **d** (int, optional) -- Max. distance (range) for grouping transcript start sites (Default: 100)

  - **e** (int, optional) -- Max. distance (range) allowed from free ends of terminal exons of reference transcripts when assessing exon accuracy (Default: 100)

  - **r** (str, optional) -- An optional "reference" annotation GFF file containing a set of known mRNAs to use as a reference for assessing the accuracy of mRNAs or gene models given in <input.gtf>

  - **run_id** (str, optional) -- An arbitrary name of the new run (which is a merge of all samples).
    - default value: magic

  - **s** (str, optional) -- Can be a multi-fasta file with all the genomic sequences or a directory containing multiple single-fasta files (one file per contig); lower case bases will be used to classify input transcripts as repeats. NOTE that must contain one fasta file per reference chromosome, and each file must be named after the chromosome, and have a .fa or .fasta extension.


**Required tools:** cuffcompare

**CPU Cores:** 1

.. index:: cufflinks

cufflinks
=========



    CuffLinks is part of the 'Cufflinks suite of tools' for
    differential expr. analysis of RNA-Seq data and their
    visualisation. This step applies the cufflinks tool which
    assembles transcriptomes from RNA-Seq data and quantifies their
    expression and produces .gtf files with these annotations.
    For details on cufflinks we refer to the author's webpage:

    http://cole-trapnell-lab.github.io/cufflinks/


**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/skipped**
  - **out/features**
  - **out/isoforms_fpkm**
  - **out/genes-fpkm**
  - **out/log_stderr**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      cufflinks [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> cufflinks;
      out_1 [label="features"];
      cufflinks -> out_1;
      out_2 [label="genes-fpkm"];
      cufflinks -> out_2;
      out_3 [label="isoforms_fpkm"];
      cufflinks -> out_3;
      out_4 [label="log_stderr"];
      cufflinks -> out_4;
      out_5 [label="skipped"];
      cufflinks -> out_5;
   }

**Options:**
  - **3-overhang-tolerance** (int, optional) -- The number of bp allowed to overhang the 3 prime end of a reference transcript when determining if an assembled transcript should be merged with it (i.e., the assembled transcript is not novel). Default: 600

  - **GTF** (bool, optional) -- Quantitate against reference transcript annotations. Use with either RABT or ab initio assembly is not supported.

  - **GTF-guide** (bool, optional) -- Tells Cufflinks to use the supplied reference annotation a GFF file to guide RABT assembly. Reference transcripts will be tiled with faux-reads to provide additional information in assembly. Output will include all reference transcripts as well as any novel genes and isoforms that are assembled.

  - **compatible-hits-norm** (bool, optional) -- With this option, Cufflinks counts only those fragments compatible with some reference transcript towards the number of mapped hits used in the FPKM denominator. This option can be combined with -N/--upper-quartile-norm. Default: FALSE

  - **frag-bias-correct** (str, optional) -- Providing Cufflinks with a multifasta file via this option instructs it to run our new bias detection and correction algorithm which can significantly improve accuracy of transcript abundance estimates. Default: NULL

  - **frag-len-mean** (int, optional) -- This is the expected (mean) fragment length. The default is 200bp. Note: Cufflinks now learns the fragment length mean for each SAM file, so using this option is no longer recommended with paired-end reads. Default: 200

  - **frag-len-std-dev** (int, optional) -- The standard deviation for the distribution on fragment lengths. The default is 80bp. Note: Cufflinks now learns the fragment length standard deviation for each SAM file, so using this option is no longer recommended with paired-end reads. Default: 80

  - **intron-overhang-tolerance** (int, optional) -- The number of bp allowed to enter the intron of a reference transcript when determining if an assembled transcript should be merged with it (i.e., the assembled transcript is not novel). Default: 50

  - **junc-alpha** (float, optional) -- The alpha value for the binomial test used during false positive spliced alignment filtration. Default: 0.001

  - **label** (str, optional) -- Cufflinks will report transfrags in GTF format, with a prefix given by this option. Default: CUFF

  - **library-norm-method** (str, optional) -- You can control how library sizes (i.e. sequencing depths) are normalized in Cufflinks and Cuffdiff. Cuffdiff has several methods that require multiple libraries in order to work. Library normalization methods supported by Cufflinks work on one library at a time. Normalization Method supported by Cufflinks: classic-fpkm (Library size factor is set to 1 - no scaling applied to FPKM values or fragment counts. Default: classic-fpkm
    - possible values: 'classic-fpkm'


  - **library-type** (str, required) -- In cases where Cufflinks cannot determine the platform and protocol used to generate input reads, you can supply this information manually, which will allow Cufflinks to infer source strand information with certain protocols. The available options are listed below. For paired-end data, we currently only support protocols where reads are point towards each other.Library type: fr-unstranded (default); examples: Standard Illumina; description: Reads from the left-most end of the fragment (in transcript coordinates) map to the transcript strand, and the right-most end maps to the opposite strand.Library type: fr-firststrand; examples: dUTP, NSR, NNSR; description: same as fr-unstranded except we enforce the rule that the right-most end of the fragment (in transcript coordinates) is the first sequenced (or only sequenced for single-end reads). Equivalently, it is assumed that only the strand generated during first strand synthesis is sequenced.Library type: fr-secondstrand; examples: Directional Illumina (Ligation), Standard SOLiD; same as fr-unstranded except we enforce the rule that the left-most end of the fragment (in transcript coordinates) is the first sequenced (or only sequenced for single-end reads). Equivalently, it is assumed that only the strand generated during second strand synthesis is sequenced. Default: fr-unstranded
    - possible values: 'ff-firststrand', 'fr-firststrand', 'ff-secondstrand', 'fr-secondstrand', 'ff-unstranded', 'fr-unstranded', 'transfrags'


  - **mask-file** (str, optional) -- Tells Cufflinks to ignore all reads that could have come from transcripts in this GTF file. We recommend including any annotated rRNA, mitochondrial transcripts other abundant transcripts you wish to ignore in your analysis in this file. Due to variable efficiency of mRNA enrichment methods and rRNA depletion kits, masking these transcripts often improves the overall robustness of transcript abundance estimates.

  - **max-bundle-frags** (int, optional) -- Sets the maximum number of fragments a locus may have before being skipped. Skipped loci are listed in skipped.gtf. Default: 500000

  - **max-bundle-length** (int, optional) -- Maximum genomic length allowed for a given bundle. Default: 3500000

  - **max-frag-multihits** (str, optional) -- Maximum number of alignments allowed per fragment. Default: unlim

  - **max-intron-length** (int, optional) -- The maximum intron length. Cufflinks will not report transcripts with introns longer than this, and will ignore SAM alignments with REF_SKIP CIGAR operations longer than this. Default: 300000

  - **max-mle-iterations** (int, optional) -- Sets the number of iterations allowed during maximum likelihood estimation of abundances. Default: 5000

  - **max-multiread-fraction** (float, optional) -- The fraction a transfrags supporting reads that may be multiply mapped to the genome. A transcript composed of more than this fraction will not be reported by the assembler. Default: 0.75 (75% multireads or more is suppressed).

  - **min-frags-per-transfrag** (int, optional) -- Assembled transfrags supported by fewer than this many aligned RNA-Seq fragments are not reported. Default: 10

  - **min-intron-length** (int, optional) -- Minimum intron size allowed in genome. Default: 50

  - **min-isoform-fraction** (float, optional) -- After calculating isoform abundance for a gene, Cufflinks filters out transcripts that it believes are very low abundance, because isoforms expressed at extremely low levels often cannot reliably be assembled, and may even be artifacts of incompletely spliced precursors of processed transcripts. This parameter is also used to filter out introns that have far fewer spliced alignments supporting them. The default is 0.1, or 10% of the most abundant isoform (the major isoform) of the gene. Range: 0.0-1.0. Default: 0.10

  - **multi-read-correct** (bool, optional) -- Tells Cufflinks to do an initial estimation procedure to more accurately weight reads mapping to multiple locations in the genome. Default: FALSE

  - **no-effective-length-correction** (bool, optional) -- Cufflinks will not employ its "effective" length normalization to transcript FPKM. Default: FALSE

  - **no-faux-reads** (bool, optional) -- This option disables tiling of the reference transcripts with faux reads. Use this if you only want to use sequencing reads in assembly but do not want to output assembled transcripts that lay within reference transcripts. All reference transcripts in the input annotation will also be included in the output. Default: FALSE

  - **no-length-correction** (bool, optional) -- Cufflinks will not normalize fragment counts by transcript length at all. Use this option when fragment count is independent of the size of the features being quantified (e.g. for small RNA libraries, where no fragmentation takes place, or 3 prime end sequencing, where sampled RNA fragments are all essentially the same length). Experimental option, use with caution. Default: FALSE

  - **no-update-check** (bool, optional) -- Turns off the automatic routine that contacts the Cufflinks server to check for a more recent version. Default: FALSE

  - **num-frag-assign-draws** (int, optional) -- Number of fragment assignment samples per generation. Default: 50

  - **num-frag-count-draws** (int, optional) -- Number of fragment generation samples. Default: 100

  - **num-threads** (int, optional) -- Number of threads used during analysis. Default: 1

  - **overhang-tolerance** (int, optional) -- The number of bp allowed to enter the intron of a transcript when determining if a read or another transcript is mappable to/compatible with it. The default is 8 bp based on the default bowtie/TopHat parameters. Default: 8

  - **overlap-radius** (int, optional) -- Transfrags that are separated by less than this distance (in bp) get merged together, and the gap is filled. Default: 50

  - **pre-mrna-fraction** (float, optional) -- Some RNA-Seq protocols produce a significant amount of reads that originate from incompletely spliced transcripts, and these reads can confound the assembly of fully spliced mRNAs. Cufflinks uses this parameter to filter out alignments that lie within the intronic intervals implied by the spliced alignments. The minimum depth of coverage in the intronic region covered by the alignment is divided by the number of spliced reads, and if the result is lower than this parameter value, the intronic alignments are ignored. The default is 15%. Range: 0.0-1.0. Default: 0.15

  - **seed** (int, optional) -- Value of random number generator seed. Default: 0

  - **small-anchor-fraction** (float, optional) -- Spliced reads with less than this percent of their length on each side of the junction are considered suspicious and are candidates for filtering prior to assembly. Default: 0.09

  - **total-hits-norm** (bool, optional) -- With this option, Cufflinks counts all fragments, including those not compatible with any reference transcript, towards the number of mapped hits used in the FPKM denominator. Default: TRUE

  - **trim-3-avgcov-thresh** (int, optional) -- Minimum average coverage required to attempt 3 prime trimming. Default: 10

  - **trim-3-dropoff-frac** (float, optional) -- The fraction of average coverage below which to trim the 3 prime end of an assembled transcript. Default: 0.1

  - **upper-quartile-norm** (bool, optional) -- DEPRECATED! Use --library-norm-method With this option, Cufflinks normalizes by the upper quartile of the number of fragments mapping to individual loci instead of the total number of sequenced fragments. This can improve robustness of differential expression calls for less abundant genes and transcripts.

  - **verbose** (bool, optional) -- Print lots of status updates and other diagnostic information. Default: FALSE


**Required tools:** cufflinks, mkdir (coreutils), mv (coreutils)

**CPU Cores:** 6

.. index:: cuffmerge

cuffmerge
=========



    CuffMerge is part of the 'Cufflinks suite of tools' for
    differential expr. analysis of RNA-Seq data and their
    visualisation. This step applies the cuffmerge tool which merges
    several Cufflinks assemblies. For details on cuffmerge we refer to
    the author's webpage:

    http://cole-trapnell-lab.github.io/cufflinks/cuffmerge/


**Input Connection**
  - **in/features**

**Output Connection**
  - **out/assemblies**
  - **out/features**
  - **out/run_log**
  - **out/log_stderr**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      cuffmerge [style=filled, fillcolor="#fce94f"];
      in_0 [label="features"];
      in_0 -> cuffmerge;
      out_1 [label="assemblies"];
      cuffmerge -> out_1;
      out_2 [label="features"];
      cuffmerge -> out_2;
      out_3 [label="log_stderr"];
      cuffmerge -> out_3;
      out_4 [label="run_log"];
      cuffmerge -> out_4;
   }

**Options:**
  - **num-threads** (int, optional) -- Use this many threads to merge assemblies.
    - default value: 6

  - **ref-gtf** (str, optional) -- A "reference" annotation GTF. The input assemblies are merged together with the reference GTF and included in the final output.

  - **ref-sequence** (str, optional) -- This argument should point to the genomic DNA sequences for the reference. If a directory, it should contain one fasta file per contig. If a multifasta file, all contigs should be present.

  - **run_id** (str, optional) -- An arbitrary name of the new run (which is a merge of all samples).
    - default value: magic


**Required tools:** cuffmerge, mkdir (coreutils), mv (coreutils), printf (coreutils)

**CPU Cores:** 6

.. index:: cutadapt

cutadapt
========



    Cutadapt finds and removes adapter sequences, primers, poly-A tails and
    other types of unwanted sequence from your high-throughput sequencing reads.

    https://cutadapt.readthedocs.org/en/stable/

    This step wraps release: cutadpat 1.5


**Input Connection**
  - **in/second_read** (optional)
  - **in/first_read**

**Output Connection**
  - **out/log_second_read** (optional)
  - **out/second_read** (optional)
  - **out/log_first_read**
  - **out/first_read**


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
      in_1 [label="second_read", style=filled, fillcolor="#a7a7a7"];
      in_1 -> cutadapt;
      out_2 [label="first_read"];
      cutadapt -> out_2;
      out_3 [label="log_first_read"];
      cutadapt -> out_3;
      out_4 [label="log_second_read", style=filled, fillcolor="#a7a7a7"];
      cutadapt -> out_4;
      out_5 [label="second_read", style=filled, fillcolor="#a7a7a7"];
      cutadapt -> out_5;
   }

**Options:**
  - **adapter-R1** (str, optional) -- Adapter sequence to be clipped off of thefirst read.

  - **adapter-R2** (str, optional) -- Adapter sequence to be clipped off of thesecond read

  - **adapter-file** (str, optional) -- File containing adapter sequences to be clipped off of the reads.

  - **adapter-type** (str, optional) -- The type of the adapter that has been used for sequencing. a: adapter ligated to the 3' end; b: adapter ligated to the 3' or 5' end (If the adapter is found within the read or overlapping the 3' end of the read, the behavior is the same as for the -a value. If the adapter overlaps the 5' end (beginning of the read), the initial portion of the read matching the adapter is trimmed, but anything that follows is kept.); g: adapter ligated to the 5' end (If the adapter sequence starts with the character '^', the adapter is 'anchored'. An anchored adapter must appear in its entirety at the 5' end of the read (it is a prefix of the read). A non-anchored adapter may appear partially at the 5' end, or it may occur within the read. If it is found within a read, the sequence preceding the adapter is also trimmed. In all cases, the adapter itself is trimmed).
    - default value: -a
    - possible values: '-a', '-g', '-b'


  - **bwa** (bool, optional) -- BWA-compatible color space output. This enables colorspace, double-encode, trim-primer, strip-f3 and suffix:'/1'.

  - **colospace** (bool, optional) -- Colorspace mode: Also trim the color that is adjacent to the found adapter.

  - **cut** (int, optional) -- Remove bases from the beginning or end of each read. If LENGTH is positive, the bases are removed from the beginning of each read. If LENGTH is negative, the bases are removed from the end of each read.

  - **dd-blocksize** (str, optional)    - default value: 2M

  - **discard-trimmed** (bool, optional) -- Discard reads that contain the adapter instead of trimming them. Also use -O in order to avoid throwing away too many randomly matching reads!

  - **discard-untrimmed** (bool, optional) -- Discard reads that do not contain the adapter.

  - **double-encode** (bool, optional) -- When in color space, double-encode colors (map 0,1,2,3,4 to A,C,G,T,N).

  - **error-rate** (float, optional) -- Maximum allowed error rate (no. of errors divided by the length of the matching region) (default: 0.1)

  - **fix_qnames** (bool, required) -- If set to true, only the leftmost string without spaces of the QNAME field of the FASTQ data is kept. This might be necessary for downstream analysis.

  - **length-tag** (str, optional) -- Search for TAG followed by a decimal number in the name of the read (description/comment field of the FASTA or FASTQ file). Replace the decimal number with the correct length of the trimmed read. For example, use --length-tag 'length=' to correct fields like 'length=123'.

  - **maq** (bool, optional) -- MAQ-compatible color space output. This enables colorspace, double-encode, trim-primer, strip-f3 and suffix:'/1'.

  - **mask-adapter** (bool, optional) -- Mask with 'N' adapter bases instead of trim (default: False)

  - **match-read-wildcards** (bool, optional) -- Allow 'N's in the read as matches to the adapter (default: False).

  - **maximum-length** (int, optional) -- Discard trimmed reads that are longer than LENGTH. Reads that are too long even before adapter removal are also discarded. In colorspace, an initial primer is not counted (default: no limit).

  - **minimum-length** (int, optional) -- Discard trimmed reads that are shorter than LENGTH. Reads that are too short even before adapter removal are also discarded. In colorspace, an initial primer is not counted (default: 0).

  - **no-indels** (bool, optional) -- Do not allow indels in the alignments, that is, allow only mismatches. This option is currently only supported for anchored 5' adapters (adapter-type: "-g" and adapter-R[1|2]: "^ADAPTER") (default: both mismatches and indels are allowed)

  - **no-trim** (bool, optional) -- Match and redirect reads to output/untrimmed-output as usual, but don't remove the adapters. (Default: False)

  - **no-zero-cap** (bool, optional) -- Do not change negative quality values to zero. Colorspace quality values of -1 would appear as spaces in the output FASTQ file. Since many tools have problems with that, negative qualities are converted to zero when trimming colorspace data. Use this option to keep negative qualities.

  - **overlap** (int, optional) -- Minimum overlap length. If the overlap between the read and the adapter is shorter than LENGTH, the read is not modified. This reduces the no. of bases trimmed purely due to short random adapter matches (default: 3).

  - **pigz-blocksize** (str, optional)    - default value: 2048

  - **prefix** (str, optional) -- Add this prefix to read names

  - **quality-base** (int, optional) -- Assume that quality values are encoded as ascii (quality + QUALITY_BASE). The default (33) is usually correct, except for reads produced by some versions of the Illumina pipeline, where this should be set to 64. (Default: 33)
    - possible values: '33', '64'


  - **quality-cutoff** (int, optional) -- Trim low-quality ends from reads before adapter removal. The algorithm is the same as the one used by BWA (Subtract CUTOFF from all qualities; compute partial sums from all indices to the end of the sequence; cut sequence at the index at which the sum is minimal) (default: 0)

  - **strip-f3** (bool, optional) -- For color space: Strip the _F3 suffix of read names

  - **strip-suffix** (str, optional) -- Remove this suffix from read names if present. Can be given multiple times.

  - **suffix** (str, optional) -- Add this suffix to read names

  - **times** (int, optional) -- Try to remove adapters at most COUNT times. Useful when an adapter gets appended multiple times (default: 1).

  - **trim-primer** (bool, optional) -- When in color space, trim primer base and the first color (which is the transition to the first nucleotide)

  - **use_reverse_complement** (bool, required) -- The reverse complement of adapter sequences 'adapter-R1' and 'adapter-R2' are used for adapter clipping.

  - **zero-cap** (bool, optional) -- Change negative quality values to zero. This is enabled by default when -c/--colorspace is also enabled. Use the above option to disable it.


**Required tools:** cat (coreutils), cutadapt, dd (coreutils), mkfifo (coreutils), pigz

**CPU Cores:** 4

.. index:: deepTools_bamCompare

deepTools_bamCompare
====================



    This tool compares two BAM files based on the number of mapped reads. To
    compare the BAM files, the genome is partitioned into bins of equal size,
    then the number of reads found in each bin is counted per file, and finally
    a summary value is reported. This value can be the ratio of the number of
    reads per bin, the log2 of the ratio, or the difference. This tool can
    normalize the number of reads in each BAM file using the SES method proposed
    by Diaz et al. (2012) "Normalization, bias correction, and peak calling for
    ChIP-seq". Statistical Applications in Genetics and Molecular Biology, 11(3).
    Normalization based on read counts is also available. The output is either a
    bedgraph or bigWig file containing the bin location and the resulting
    comparison value. By default, if reads are paired, the fragment length
    reported in the BAM file is used. Each mate, however, is treated
    independently to avoid a bias when a mixture of concordant and discordant
    pairs is present. This means that each end will be extended to match the
    fragment length.

    http://deeptools.readthedocs.io/en/latest/content/tools/bamCompare.html

    Usage example::

        bamCompare -b1 treatment.bam -b2 control.bam -o log2ratio.bw


**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/ucsc-tracks**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      deepTools_bamCompare [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> deepTools_bamCompare;
      out_1 [label="ucsc-tracks"];
      deepTools_bamCompare -> out_1;
   }

**Options:**
  - **binSize** (int, optional) -- Size of the bins, in bases, for the output of the bigwig/bedgraph file. (default: 50)

  - **blackListFileName** (str, optional) -- A BED or GTF file containing regions that should be excluded from all analyses. Currently this works by rejecting genomic chunks that happen to overlap an entry. Consequently, for BAM files, if a read partially overlaps a blacklisted region or a fragment spans over it, then the read/fragment might still be considered. Please note that you should adjust the effective genome size, if relevant. (default: None)

  - **centerReads** (bool, optional) -- By adding this option, reads are centered with respect to the fragment length. For paired-end data, the read is centered at the fragment length defined by the two ends of the fragment. For single-end data, the given fragment length is used. This option is useful to get a sharper signal around enriched regions. (default: False)

  - **extendReads** (int, optional) -- This parameter allows the extension of reads to fragment size. If set, each read is extended, without exception. *NOTE*: This feature is generally NOT recommended for spliced-read data, such as RNA-seq, as it would extend reads over skipped regions. *Single-end*: Requires a user specified value for the final fragment length. Reads that already exceed this fragment length will not be extended. *Paired-end*: Reads with mates are always extended to match the fragment size defined by the two read mates. Unmated reads, mate reads that map too far apart (>4x fragment length) or even map to different chromosomes are treated like single-end reads. The input of a fragment length value is optional. If no value is specified, it is estimated from the data (mean of the fragment size of all mate reads). (default: False)

  - **ignoreDuplicates** (bool, optional) -- If set, reads that have the same orientation and start position will be considered only once. If reads are paired, the mate's position also has to coincide to ignore a read. (default: False)

  - **ignoreForNormalization** (list, optional) -- A list of space-delimited chromosome names containing those chromosomes that should be excluded for computing the normalization. This is useful when considering samples with unequal coverage across chromosomes, like male samples. An usage examples is --ignoreForNormalization chrX chrM. (default: None)

  - **maxFragmentLength** (int, optional) -- The maximum fragment length needed for read/pair inclusion. A value of 0 disables filtering and is needed for including single-end and orphan reads. (default: 0)

  - **minFragmentLength** (int, optional) -- The minimum fragment length needed for read/pair inclusion. Note that a value other than 0 will exclude all single-end reads. This option is primarily useful in ATACseq experiments, for filtering mono- or di-nucleosome fragments. (default: 0)

  - **minMappingQuality** (int, optional) -- If set, only reads that have a mapping quality score of at least this are considered. (default: None)

  - **normalizeTo1x** (int, optional) -- Report read coverage normalized to 1x sequencing depth (also known as Reads Per Genomic Content (RPGC)). Sequencing depth is defined as: (total number of mapped reads * fragment length) / effective genome size. The scaling factor used is the inverse of the sequencing depth computed for the sample to match the 1x coverage. To use this option, the effective genome size has to be indicated after the option. The effective genome size is the portion of the genome that is mappable. Large fractions of the genome are stretches of NNNN that should be discarded. Also, if repetitive regions were not included in the mapping of reads, the effective genome size needs to be adjusted accordingly. Common values are: mm9:2,150,570,000; hg19:2,451,960,000; dm3:121,400,000 and ce10:93,260,000. See Table 2 of http://www.plosone.org/article/info:doi/10.1371/journal.pone.0030377 or http://www.nature.com/nbt/journal/v27/n1/fig_tab/nbt.1518_T1.html for several effective genome sizes. (default: None)

  - **normalizeUsingRPKM** (bool, optional) -- Use Reads Per Kilobase per Million reads to normalize the number of reads per bin. The formula is: RPKM (per bin) = number of reads per bin / ( number of mapped reads (in millions) * bin length (kb) ). Each read is considered independently,if you want to only count either of the mate pairs in paired-end data, use the --samFlag option. (default: False)

  - **numberOfSamples** (int, optional) -- *Only relevant when SES is chosen for the scaleFactorsMethod.* Number of samplings taken from the genome to compute the scaling factors. (default: 100000.0)

  - **outFileFormat** (str, required) -- Output file type. Either "bigwig" or "bedgraph". (default: "bigwig")
    - possible values: 'bigwig', 'bedgraph'


  - **pseudocount** (float, optional) -- small number to avoid x/0. Only useful together with --ratio log2 or --ratio ratio . (default: 1)

  - **ratio** (str, optional) -- The default is to output the log2ratio of the two samples. The reciprocal ratio returns the negative of the inverse of the ratio if the ratio is less than 0. The resulting values are interpreted as negative fold changes. *NOTE*: Only with --ratio subtract can --normalizeTo1x or --normalizeUsingRPKM be used. Instead of performing a computation using both files, the scaled signal can alternatively be output for the first or second file using the '--ratio first' or '--ratio second' (default: log2)
    - possible values: 'log2', 'ratio', 'subtract', 'add', 'mean', 'reciprocal_ratio', 'first,second'


  - **region** (str, optional) -- Region of the genome to limit the operation to - this is useful when testing parameters to reduce the computing time. The format is chr:start:end, for example --region chr10 or --region chr10:456700:891000. (default: None)

  - **samFlagExclude** (int, optional) -- Exclude reads based on the SAM flag. For example, to get only reads that map to the forward strand, use --samFlagExclude 16, where 16 is the SAM flag for reads that map to the reverse strand. (default: None)

  - **samFlagInclude** (int, optional) -- Include reads based on the SAM flag. For example, to get only reads that are the first mate, use a flag of 64. This is useful to count properly paired reads only once, as otherwise the second mate will be also considered for the coverage. (default: None)

  - **sampleLength** (int, optional) -- *Only relevant when SES is chosen for the scaleFactorsMethod.* To compute the SES, specify the length (in bases) of the regions (see --numberOfSamples) that will be randomly sampled to calculate the scaling factors. If you do not have a good sequencing depth for your samples consider increasing the sampling regions' size to minimize the probability that zero-coverage regions are used. (default: 1000)

  - **samples** (list, required) -- List of lists with two elements. Each element has to be the name of a run. Each run has to provide a SINGLE BAM file. Both BAM files are compared using deepTools bamCompare command.

  - **scaleFactors** (str, optional) -- Set this parameter manually to avoid the computation of scaleFactors. The format is scaleFactor1:scaleFactor2. For example, --scaleFactor 0.7:1 will cause the first BAM file tobe multiplied by 0.7, while not scaling the second BAM file (multiplication with 1). (default: None)

  - **scaleFactorsMethod** (str, optional) -- Method to use to scale the samples. (default: readCount)
    - possible values: 'readCount', 'SES'


  - **skipNonCoveredRegions** (bool, optional) -- This parameter determines if non-covered regions (regions without overlapping reads) in a BAM file should be skipped. The default is to treat those regions as having a value of zero. The decision to skip non-covered regions depends on the interpretation of the data. Non-covered regions may represent, for example, repetitive regions that should be skipped. (default: False)

  - **smoothLength** (int, optional) -- The smooth length defines a window, larger than the binSize, to average the number of reads. For example, if the --binSize is set to 20 and the --smoothLength is set to 60, then, for each bin, the average of the bin and its left and right neighbors is considered. Any value smaller than --binSize will be ignored and no smoothing will be applied. (default: None)


**Required tools:** bamCompare

**CPU Cores:** 10

.. index:: deepTools_bamPEFragmentSize

deepTools_bamPEFragmentSize
===========================



    bamPEFragmentSize tool calculates the fragment sizes for read pairs given a
    BAM file from paired-end sequencing.Several regions are sampled depending on
    the size of the genome and number of processors to estimate thesummary
    statistics on the fragment lengths. Properly paired reads are preferred for
    computation, i.e., it will only use discordant pairs if no concordant
    alignments overlap with a given region. The default setting simply prints
    the summary statistics to the screen.

    http://deeptools.readthedocs.io/en/latest/content/tools/bamPEFragmentSize.html

    Usage example::

        bamPEFragmentSize [-h] [--bamfiles bam files [bam files ...]]
                          [--histogram FILE] [--numberOfProcessors INT]
                          [--samplesLabel SAMPLESLABEL [SAMPLESLABEL ...]]
                          [--plotTitle PLOTTITLE]
                          [--maxFragmentLength MAXFRAGMENTLENGTH] [--logScale]
                          [--binSize INT] [--distanceBetweenBins INT]
                          [--blackListFileName BED file] [--verbose] [--version]


**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/fragment_size_stats**
  - **out/fragment_size_plots**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      deepTools_bamPEFragmentSize [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> deepTools_bamPEFragmentSize;
      out_1 [label="fragment_size_plots"];
      deepTools_bamPEFragmentSize -> out_1;
      out_2 [label="fragment_size_stats"];
      deepTools_bamPEFragmentSize -> out_2;
   }

**Options:**
  - **binSize** (int, optional) -- Length in bases of the window used to sample the genome. (default 1000)

  - **blackListFileName** (str, optional) -- A BED file containing regions that should be excluded from all analyses. Currently this works by rejecting genomic chunks that happen to overlap an entry. Consequently, for BAM files, if a read partially overlaps a blacklisted region or a fragment spans over it, then the read/fragment might still be considered. Please note that you should adjust the effective genome size, if relevant. (default: None)

  - **distanceBetweenBins** (int, optional) -- To reduce the computation time, not every possible genomic bin is sampled. This option allows you to set the distance between bins actually sampled from. Larger numbers are sufficient for high coverage samples, while smaller values are useful for lower coverage samples. Note that if you specify a value that results in too few (<1000) reads sampled, the value will be decreased. (default 1000000)

  - **histogram** (bool, optional) -- If set saves a .png file with a histogram of fragment length distribution for each run.

  - **logScale** (bool, optional) -- Plot on the log scale

  - **maxFragmentLength** (int, optional) -- The maximum fragment length in the histogram. A value of 0 (the default) indicates to use twice the mean fragment length

  - **samples** (dict, optional) -- Dictionary with IDs of new runs as keys and lists of sample names as values. For each sample name a BAM file is expected to be the input from upstream steps. If not provided this step calculates summary statistics for each input file.


**Required tools:** bamPEFragmentSize

**CPU Cores:** 10

.. index:: deepTools_multiBamSummary

deepTools_multiBamSummary
=========================



    This step computes the read coverages for genomic regions for every BAM
    input file using multiBamSummary. For downstream analysis such as
    'plotCorrelation' or 'plotPCA' you need to merge the output files.

    The analysis can be performed for the
    entire genome by running the program in 'bins' mode. If you want to count
    the read coverage for specific regions only, use the BED-file mode instead.
    The standard output of multiBamSummary is a compressed numpy array (.npz).
    It can be directly used to calculate and visualize pairwise correlation
    values between the read coverages using the tool 'plotCorrelation'.
    Similarly, plotPCA can be used for principal component analysis of the read
    coverages using the .npz file. Note that using a single bigWig file is only
    recommended if you want to produce a bedGraph file (i.e., with the
    --outRawCounts option; the default output file cannot be used by ANY
    deepTools program if only a single file was supplied!).

    http://deeptools.readthedocs.io/en/latest/content/tools/multiBamSummary.html

    Usage example::

        multiBamSummary [-h] [--version]  ...


**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/read-coverage**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      deepTools_multiBamSummary [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> deepTools_multiBamSummary;
      out_1 [label="read-coverage"];
      deepTools_multiBamSummary -> out_1;
   }

**Options:**
  - **bed-file** (list, optional) -- BED file that contains all regions that should be considered for the coverage analysis. If this option is set "multiBamSummary" is executed with "BED-file" subcommand, otherwise with "bins" subcommand.

  - **binSize** (int, optional) -- Length in bases of the window used to sample the genome. (default: 10000)

  - **blackListFileName** (str, optional) -- A BED or GTF file containing regions that should be excluded from all analyses. Currently this works by rejecting genomic chunks that happen to overlap an entry. Consequently, for BAM files, if a read partially overlaps a blacklisted region or a fragment spans over it, then the read/fragment might still be considered. Please note that you should adjust the effective genome size, if relevant. (default: None)

  - **centerReads** (bool, optional) -- By adding this option, reads are centered with respect to the fragment length. For paired-end data, the read is centered at the fragment length defined by the two ends of the fragment. For single-end data, the given fragment length is used. This option is useful to get a sharper signal around enriched regions. (default: False)

  - **distanceBetweenBins** (int, optional) -- By default, multiBamSummary considers consecutive bins of the specified --binSize. However, to reduce the computation time, a larger distance between bins can by given. Larger distances result in fewer bins considered. (default: 0)

  - **exonID** (str, optional) -- When a GTF file is used to provide regions, only entries with this value as their feature (column 2) will be processed as exons. CDS would be another common value for this. (default: exon)

  - **extendReads** (bool, optional) -- This parameter allows the extension of reads to fragment size. If set, each read is extended, without exception. *NOTE*: This feature is generally NOT recommended for spliced-read data, such as RNA-seq, as it would extend reads over skipped regions. *Single-end*: Requires a user specified value for the final fragment length. Reads that already exceed this fragment length will not be extended. *Paired-end*: Reads with mates are always extended to match the fragment size defined by the two read mates. Unmated reads, mate reads that map too far apart (>4x fragment length) or even map to different chromosomes are treated like single-end reads. The input of a fragment length value is optional. If no value is specified, it is estimated from the data (mean of the fragment size of all mate reads). (default: False)

  - **ignoreDuplicates** (bool, optional) -- If set, reads that have the same orientation and start position will be considered only once. If reads are paired, the mate's position also has to coincide to ignore a read. (default: False)

  - **maxFragmentLength** (int, optional) -- The maximum fragment length needed for read/pair inclusion. A value of 0 disables filtering and is needed for including single-end and orphan reads. (default: 0)

  - **metagene** (bool, optional) -- When either a BED12 or GTF file are used to provide regions, perform the computation on the merged exons, rather than using the genomic interval defined by the 5-prime and 3-prime most transcript bound (i.e., columns 2 and 3 of a BED file). If a BED3 or BED6 file is used as input, then columns 2 and 3 are used as an exon. (default: False)

  - **minFragmentLength** (int, optional) -- The minimum fragment length needed for read/pair inclusion. Note that a value other than 0 will exclude all single-end reads. This option is primarily useful in ATACseq experiments, for filtering mono- or di-nucleosome fragments. (default: 0)

  - **minMappingQuality** (int, optional) -- If set, only reads that have a mapping quality score of at least this are considered. (default: None)

  - **outRawCounts** (bool, optional) -- Save the counts per region to a tab-delimited file. (default: False)

  - **region** (str, optional) -- Region of the genome to limit the operation to - this is useful when testing parameters to reduce the computing time. The format is chr:start:end, for example --region chr10 or --region chr10:456700:891000. (default: None)

  - **samFlagExclude** (int, optional) -- Exclude reads based on the SAM flag. For example, to get only reads that map to the forward strand, use --samFlagExclude 16, where 16 is the SAM flag for reads that map to the reverse strand. (default: None)

  - **samFlagInclude** (int, optional) -- Include reads based on the SAM flag. For example, to get only reads that are the first mate, use a flag of 64. This is useful to count properly paired reads only once, as otherwise the second mate will be also considered for the coverage. (default: None)

  - **transcriptID** (str, optional) -- When a GTF file is used to provide regions, only entries with this value as their feature (column 2) will be processed as transcripts. (default: transcript)

  - **transcript_id_designator** (str, optional) -- Each region has an ID (e.g., ACTB) assigned to it, which for BED files is either column 4 (if it exists) or the interval bounds. For GTF files this is instead stored in the last column as a key:value pair (e.g., as 'transcript_id "ACTB"', for a key of transcript_id and a value of ACTB). In some cases it can be convenient to use a different identifier. To do so, set this to the desired key. (default: transcript_id)


**Required tools:** multiBamSummary

**CPU Cores:** 10

.. index:: deepTools_plotFingerprint

deepTools_plotFingerprint
=========================



    This tool samples indexed BAM files and plots a profile of cumulative
    read coverages for each. All reads overlapping a window (bin) of the
    specified length are counted; these counts are sorted and the cumulative
    sum is finally plotted.

    Usage example::

       plotFingerprint -b treatment.bam control.bam -plot fingerprint.png


**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/counts**
  - **out/plots**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      deepTools_plotFingerprint [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> deepTools_plotFingerprint;
      out_1 [label="counts"];
      deepTools_plotFingerprint -> out_1;
      out_2 [label="plots"];
      deepTools_plotFingerprint -> out_2;
   }

**Options:**
  - **JSDsample** (str, optional) -- Reference sample against which to compute the Jensen-Shannon distance and the CHANCE statistics. If this is not specified, then these will not be calculated. If --outQualityMetrics is not specified then this will be ignored. The Jensen-Shannon implementation is based on code from Sitanshu Gakkhar at BCGSC. The CHANCE implementation is based on code from Matthias Haimel. (default: None)

  - **binSize** (int, optional) -- Window size in base pairs to sample the genome. (default: 500)

  - **blackListFileName** (str, optional) -- A BED or GTF file containing regions that should be excluded from all analyses. Currently this works by rejecting genomic chunks that happen to overlap an entry. Consequently, for BAM files, if a read partially overlaps a blacklisted region or a fragment spans over it, then the read/fragment might still be considered. Please note that you should adjust the effective genome size, if relevant. (default: None)

  - **centerReads** (bool, optional) -- By adding this option, reads are centered with respect to the fragment length. For paired-end data, the read is centered at the fragment length defined by the two ends of the fragment. For single-end data, the given fragment length is used. This option is useful to get a sharper signal around enriched regions. (default: False)

  - **extendReads** (int, optional) -- This parameter allows the extension of reads to fragment size. If set, each read is extended, without exception. *NOTE*: This feature is generally NOT recommended for spliced-read data, such as RNA-seq, as it would extend reads over skipped regions. *Single-end*: Requires a user specified value for the final fragment length. Reads that already exceed this fragment length will not be extended. *Paired-end*: Reads with mates are always extended to match the fragment size defined by the two read mates. Unmated reads, mate reads that map too far apart (>4x fragment length) or even map to different chromosomes are treated like single-end reads. The input of a fragment length value is optional. If no value is specified, it is estimated from the data (mean of the fragment size of all mate reads). (default: False)

  - **ignoreDuplicates** (bool, optional) -- If set, reads that have the same orientation and start position will be considered only once. If reads are paired, the mate's position also has to coincide to ignore a read. (default: False)

  - **maxFragmentLength** (int, optional) -- The maximum fragment length needed for read/pair inclusion. A value of 0 disables filtering and is needed for including single-end and orphan reads. (default: 0)

  - **minFragmentLength** (int, optional) -- The minimum fragment length needed for read/pair inclusion. Note that a value other than 0 will exclude all single-end reads. This option is primarily useful in ATACseq experiments, for filtering mono- or di-nucleosome fragments. (default: 0)

  - **minMappingQuality** (int, optional) -- If set, only reads that have a mapping quality score of at least this are considered. (default: None)

  - **numberOfSamples** (int, optional) -- Number of bins that sampled from the genome, for which the overlapping number of reads is computed. (default: 500000.0)

  - **outQualityMetrics** (str, optional) -- Quality metrics can optionally be output to this file. The file will have one row per input BAM file and columns containing a number of metrics. Please see the online documentation for a longer explanation: http://deeptools.readthedocs.io/en/latest/content/feature/plotFingerprint_QC_metrics.html. (default: None)

  - **plotFileFormat** (str, required) -- File ending of the output figure. It will be used to determine the image format.
    - possible values: 'png', 'eps', 'pdf', 'svg'


  - **region** (str, optional) -- Region of the genome to limit the operation to - this is useful when testing parameters to reduce the computing time. The format is chr:start:end, for example --region chr10 or --region chr10:456700:891000. (default: None)

  - **samFlagExclude** (int, optional) -- Exclude reads based on the SAM flag. For example, to get only reads that map to the forward strand, use --samFlagExclude 16, where 16 is the SAM flag for reads that map to the reverse strand. (default: None)

  - **samFlagInclude** (int, optional) -- Include reads based on the SAM flag. For example, to get only reads that are the first mate, use a flag of 64. This is useful to count properly paired reads only once, as otherwise the second mate will be also considered for the coverage. (default: None)

  - **samples** (list, required) -- List of lists with run names. Each element has to be the name of a run. Each run has to provide a SINGLE BAM file. All BAM files are plotted and counted using deepTools plotFingerprint command.

  - **skipZeros** (bool, optional) -- If set, then regions with zero overlapping reads for *all* given BAM files are ignored. This will result in a reduced number of read counts than that specified in --numberOfSamples (default: False)


**Required tools:** plotFingerprint

**CPU Cores:** 10

.. index:: discardLargeSplitsAndPairs

discardLargeSplitsAndPairs
==========================



    discardLargeSplitsAndPairs reads SAM formatted alignments of the
    mapped reads. It discards all split reads that skip more than
    splits_N nucleotides in their alignment to the ref genome. In
    addition, all read pairs that are mapped to distant region such
    that the final template will exceed N_mates nucleotides will also
    be discarded. All remaining reads are returned in SAM format. The
    discarded reads are also collected in a SAM formatted file and a
    statistic is returned.

**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/log**
  - **out/alignments**
  - **out/stats**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      discardLargeSplitsAndPairs [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> discardLargeSplitsAndPairs;
      out_1 [label="alignments"];
      discardLargeSplitsAndPairs -> out_1;
      out_2 [label="log"];
      discardLargeSplitsAndPairs -> out_2;
      out_3 [label="stats"];
      discardLargeSplitsAndPairs -> out_3;
   }

**Options:**
  - **M_mates** (str, required) -- Size of template (in nucleotides) that would arise from a read pair. Read pairs that exceed this value are discarded. 

  - **N_splits** (str, required) -- Size of the skipped region within a split read (in nucleotides). Split Reads that skip more nt than this value are discarded.


**Required tools:** dd (coreutils), pigz, samtools

**CPU Cores:** 4

.. index:: fastq_screen

fastq_screen
============



    Fastq Screen (ver. 0.11.*)

    FastQ Screen allows you to screen a library of sequences in FastQ format
    against a set of sequence databases so you can see if the composition of
    the library matches with what you expect.


**Input Connection**
  - **in/first_read**

**Output Connection**
  - **out/fqc_image**
  - **out/tagged_filter** (optional)
  - **out/fastq_screen.conf** (optional)
  - **out/fqc_report**
  - **out/fqc_html**
  - **out/log_stdout**
  - **out/log_stderr**
  - **out/tagged** (optional)


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      fastq_screen [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> fastq_screen;
      out_1 [label="fastq_screen.conf", style=filled, fillcolor="#a7a7a7"];
      fastq_screen -> out_1;
      out_2 [label="fqc_html"];
      fastq_screen -> out_2;
      out_3 [label="fqc_image"];
      fastq_screen -> out_3;
      out_4 [label="fqc_report"];
      fastq_screen -> out_4;
      out_5 [label="log_stderr"];
      fastq_screen -> out_5;
      out_6 [label="log_stdout"];
      fastq_screen -> out_6;
      out_7 [label="tagged", style=filled, fillcolor="#a7a7a7"];
      fastq_screen -> out_7;
      out_8 [label="tagged_filter", style=filled, fillcolor="#a7a7a7"];
      fastq_screen -> out_8;
   }

**Options:**
  - **config** (str, optional) -- Manually specify a location for the configuration.

  - **cores** (int, required)    - default value: 10

  - **databases** (dict, optional) -- Manually specify a location for the configuration. E.g.: fastq_screen: databases: Human: /path/to/human/bowtie2/index

  - **keep config** (bool, optional) -- Keep the generated fastq_screen.conf for each run (only if databases are specified).

  - **nohits** (bool, optional) -- Writes to a file the sequences that did not map to any of the specified genomes. This option is equivalent to specifying --tag --filter 0000 (number of zeros corresponds to the number of genomes screened). By default the whole input file will be mapped, unless overridden by --subset.

  - **subset** (int, optional) -- Don't use the whole sequence file, but create a temporary dataset of this specified number of reads. The dataset created will be of approximately (within a factor of 2) of this size. If the real dataset is smaller than twice the specified size then the whole dataset will be used. Subsets will be taken evenly from throughout the whole original dataset. By Default FastQ Screen runs with this parameter set to 100,000. To process an entire dataset however, adjust --subset to 0.


**Required tools:** bowtie2, fastq_screen, mv (coreutils), printf (coreutils), rm (coreutils)

**CPU Cores:** 10

.. index:: fastqc

fastqc
======



    The fastqc step  is a wrapper for the fastqc tool. It generates some quality
    metrics for fastq files. For this specific instance only the zip archive is
    preserved.

    http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

    Tested fastqc release: 0.11.2

**Input Connection**
  - **in/second_read** (optional)
  - **in/first_read**

**Output Connection**
  - **out/first_read_fastqc_report**
  - **out/second_read_fastqc_report** (optional)
  - **out/second_read_log_stderr** (optional)
  - **out/first_read_log_stderr**
  - **out/second_read_fastqc_report_webpage** (optional)
  - **out/first_read_fastqc_report_webpage**


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
      in_1 [label="second_read", style=filled, fillcolor="#a7a7a7"];
      in_1 -> fastqc;
      out_2 [label="first_read_fastqc_report"];
      fastqc -> out_2;
      out_3 [label="first_read_fastqc_report_webpage"];
      fastqc -> out_3;
      out_4 [label="first_read_log_stderr"];
      fastqc -> out_4;
      out_5 [label="second_read_fastqc_report", style=filled, fillcolor="#a7a7a7"];
      fastqc -> out_5;
      out_6 [label="second_read_fastqc_report_webpage", style=filled, fillcolor="#a7a7a7"];
      fastqc -> out_6;
      out_7 [label="second_read_log_stderr", style=filled, fillcolor="#a7a7a7"];
      fastqc -> out_7;
   }

**Options:**
  - **adapters** (str, optional) -- Specifies a non-default file which contains the list of adapter sequences which will be explicity searched against the library. The file must contain sets of named adapters in the form name[tab]sequence. Lines prefixed with a hash will be ignored.

  - **casava** (bool, optional) -- Files come from raw casava output. Files in the same sample group (differing only by the group number) will be analysed as a set rather than individually. Sequences with the filter flag set in the header will be excluded from the analysis. Files must have the same names given to them by casava (including being gzipped and ending with .gz) otherwise they won't be grouped together correctly.

  - **contaminants** (str, optional) -- Specifies a non-default file which contains the list of contaminants to screen overrepresented sequences against. The file must contain sets of named contaminants in the form name[tab]sequence. Lines prefixed with a hash will be ignored.

  - **dd-blocksize** (str, optional)    - default value: 2M

  - **dir** (str, optional) -- Selects a directory to be used for temporary files written when generating report images. Defaults to system temp directory if not specified.

  - **format** (str, optional) -- Bypasses the normal sequence file format detection and forces the program to use the specified format. Valid formats are bam,sam, bam_mapped,sam_mapped and fastq
    - possible values: 'bam', 'sam', 'bam_mapped', 'sam_mapped', 'fastq'


  - **java** (str, optional) -- Provides the full path to the java binary you want to use to launch fastqc. If not supplied then java is assumed to be in your path.

  - **kmers** (int, optional) -- Specifies the length of Kmer to look for in the Kmer content module. Specified Kmer length must be between 2 and 10. Default length is 7 if not specified.

  - **limits** (str, optional) -- Specifies a non-default file which contains a set of criteria which will be used to determine the warn/error limits for the various modules. This file can also be used to selectively remove some modules from the output all together. The format needs to mirror the default limits.txt file found in the Configuration folder.

  - **nofilter** (bool, optional) -- If running with --casava then do not remove read flagged by casava as poor quality when performing the QC analysis.

  - **nogroup** (bool, optional) -- Disable grouping of bases for reads >50bp. All reports will show data for every base in the read. WARNING: Using this option will cause fastqc to crash and burn if you use it on really long reads, and your plots may end up a ridiculous size. You have been warned!

  - **pigz-blocksize** (str, optional)    - default value: 2048

  - **threads** (int, optional) -- Specifies the number of files which can be processed simultaneously. Each thread will be allocated 250MB of memory so you should not run more threads than your available memory will cope with, and not more than 6 threads on a 32 bit machine


**Required tools:** fastqc, mkdir (coreutils), mv (coreutils)

**CPU Cores:** 4

.. index:: fastqsample

fastqsample
===========



    wrapper class for fastq-sample
    sample random reads from a fastq file
    http://homes.cs.washington.edu/~dcjones/fastq-tools/fastq-sample.html

    for a specific seed the subsampling process will ever produce the
    same order of positions so the connections between R1 and R2 remains
    (paired end)

**Input Connection**
  - **in/second_read** (optional)
  - **in/first_read**

**Output Connection**
  - **out/second_read** (optional)
  - **out/first_read**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      fastqsample [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> fastqsample;
      in_1 [label="second_read", style=filled, fillcolor="#a7a7a7"];
      in_1 -> fastqsample;
      out_2 [label="first_read"];
      fastqsample -> out_2;
      out_3 [label="second_read", style=filled, fillcolor="#a7a7a7"];
      fastqsample -> out_3;
   }

**Options:**
  - **c** (str, optional) -- Output reads not included in the random sample to a file (or files) with the given prefix. By default, these reads are not output.

  - **n** (int, optional) -- The number of reads to sample and output
    - default value: 1000

  - **o** (str, optional) -- The filename prefix to which output should be written. If single-end data is being sampled, the output file is [PREFIX].fastq, and with paired-end, [PREFIX].1.fastq and [PREFIX].2.fastq

  - **p** (float, optional) -- The number of reads to sample in terms of the proportion of total reads. If sampling with replacement, this number may be greater than 1.0

  - **r** (str, optional) -- Sample with replacement

  - **s** (str, optional) -- Seed the random number generator. Using the same seed on the same data set will produce the same random sample.
    - default value: 1234


**Required tools:** fastq-sample, mv (coreutils), pigz, rm (coreutils)

**CPU Cores:** 1

.. index:: fastx_quality_stats

fastx_quality_stats
===================



    fastx_quality_stats generates a text file containing quality information
    of the input fastq data.

    Documentation::

        http://hannonlab.cshl.edu/fastx_toolkit/

    Tested fastqc release: 0.0.13

**Input Connection**
  - **in/second_read**
  - **in/first_read**

**Output Connection**
  - **out/second_read_quality_stats**
  - **out/first_read_quality_stats**


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
  - **dd-blocksize** (str, optional)    - default value: 2M

  - **new_output_format** (bool, optional) -- New output format (with more information per nucleotide/cycle).

  - **pigz-blocksize** (str, optional)    - default value: 2048

  - **quality** (int, optional)    - default value: 33


**Required tools:** cat (coreutils), dd (coreutils), fastx_quality_stats, mkfifo (coreutils), pigz

**CPU Cores:** 4

.. index:: fastx_reverse_complement

fastx_reverse_complement
========================



    wrapper class for fastx_reverse_complement from fastx toolkit
    creates reverse complement of fasta and fastq files.
    http://hannonlab.cshl.edu/fastx_toolkit/

**Input Connection**
  - **in/fastx**

**Output Connection**
  - **out/fastx**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      fastx_reverse_complement [style=filled, fillcolor="#fce94f"];
      in_0 [label="fastx"];
      in_0 -> fastx_reverse_complement;
      out_1 [label="fastx"];
      fastx_reverse_complement -> out_1;
   }

**Options:**
  - **prefix** (str, optional) -- Add Prefix to sample name


**Required tools:** cat (coreutils), fastx_reverse_complement, pigz

**CPU Cores:** 1

.. index:: feature_counts

feature_counts
==============



    comment here

**Input Connection**
  - **in/feature-file**
  - **in/alignments**

**Output Connection**
  - **out/counts**
  - **out/summary**
  - **out/log_stdout**
  - **out/log_stderr**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      feature_counts [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> feature_counts;
      in_1 [label="feature-file"];
      in_1 -> feature_counts;
      out_2 [label="counts"];
      feature_counts -> out_2;
      out_3 [label="log_stderr"];
      feature_counts -> out_3;
      out_4 [label="log_stdout"];
      feature_counts -> out_4;
      out_5 [label="summary"];
      feature_counts -> out_5;
   }

**Options:**
  - **A** (str, optional) -- Specify the name of a file including aliases of chromosome names. The file should be a comma delimited text file that includes two columns. The first column gives the chromosome names used in the annotation and the second column gives the chromosome names used by reads. This file should not contain header lines. Names included in this file are case sensitive.

  - **B** (bool, optional) -- If specified, only fragments that have both ends successfully aligned will be considered for summarization. This option is only applicable for paired-end reads.

  - **C** (bool, optional) -- If specified, the chimeric fragments (those fragments that have their two ends aligned to different chromosomes) will NOT be included for summarization. This option is only applicable for paired-end read data.

  - **D** (int, optional) -- Maximum fragment/template length, 600 by default.

  - **F** (str, optional) -- Specify the format of the annotation file. Acceptable formats include 'GTF' and 'SAF'. 'GTF' by default. Please refer to the users guide for SAF annotation format.

  - **M** (bool, optional) -- If specified, multi-mapping reads/fragments will be counted (ie. a multi-mapping read will be counted up to N times if it has N reported mapping locations). The program uses the 'NH' tag to find multi-mapping reads.

  - **O** (bool, optional) -- If specified, reads (or fragments if -p is specified) will be allowed to be assigned to more than one matched meta-feature (or feature if -f is specified).

  - **P** (bool, optional) -- If specified, paired-end distance will be checked when assigning fragments to meta-features or features. This option is only applicable when -p is specified. The distance thresholds should be specified using -d and -D options.

  - **Q** (int, optional) -- The minimum mapping quality score a read must satisfy in order to be counted. For paired-end reads, at least one end should satisfy this criteria. 0 by default.

  - **R** (bool, optional) -- Output read counting result for each read/fragment. For each input read file, read counting results for reads/fragments will be saved to a tab-delimited file that contains four columns including read name, status(assigned or the reason if not assigned), name of target feature/meta-feature and number of hits if the read/fragment is counted multiple times. Name of the file is the same as name of the input read file except a suffix '.featureCounts' is added.

  - **T** (int, optional) -- Number of the threads. 1 by default.

  - **a** (str, required) -- Give the name of the annotation file. The program assumes hat the provided annotation file is in GTF format. Use -F option to specify other annotation formats.

  - **cores** (int, required)    - default value: 3

  - **countSplitAlignmentsOnly** (bool, optional) -- If specified, only split alignments (CIGAR strings containing letter 'N') will be counted. All the other alignments will be ignored. An example of split alignments is the exon-spanning reads in RNA-seq data.

  - **d** (int, optional) -- Minimum fragment/template length, 50 by default.

  - **f** (bool, optional) -- If specified, read summarization will be performed at the feature level (eg. exon level). Otherwise, it is performed at meta-feature level (eg. gene level).

  - **g** (str, optional) -- Specify the attribute type used to group features (eg. exons) into meta-features (eg. genes), when GTF annotation is provided. 'gene_id' by default. This attribute type is usually the gene identifier. This argument is useful for the meta-feature level summarization.

  - **ignoreDup** (bool, optional) -- If specified, reads that were marked as duplicates will be ignored. Bit Ox400 in FLAG field of SAM/BAM file is used for identifying duplicate reads. In paired end data, the entire read pair will be ignored if at least one end is found to be a duplicate read.

  - **minReadOverlap** (int, optional) -- Specify the minimum number of overlapped bases required to assign a read to a feature. 1 by default. Negative values are permitted, indicating a gap being allowed between a read and a feature.

  - **o** (str, optional) -- Give the name of the output file. The output file contains the number of reads assigned to each meta-feature (or each feature if -f is specified). A meta-feature is the aggregation of features, grouped by using gene identifiers. Please refer to the users guide for more details.
    - default value: counts.txt

  - **p** (bool, optional) -- If specified, fragments (or templates) will be counted instead of reads. This option is only applicable for paired-end reads. The two reads from the same fragment must be adjacent to each other in the provided SAM/BAM file.

  - **primary** (bool, optional) -- If specified, only primary alignments will be counted. Primary and secondary alignments are identified using bit 0x100 in the Flag field of SAM/BAM files. All primary alignments in a dataset will be counted no matter they are from multi-mapping reads or not ('-M' is ignored).

  - **read2pos** (int, optional) -- The read is reduced to its 5' most base or 3' most base. Read summarization is then performed based on the single base which the read is reduced to. By default, no read reduction will be performed.

  - **readExtension3** (int, optional) -- Reads are extended upstream by <int> bases from their 3' end. 0 by default.

  - **readExtension5** (int, optional) -- Reads are extended upstream by <int> bases from their 5' end. 0 by default.

  - **s** (int, required) -- Indicate if strand-specific read counting should be performed. It has three possible values: 0 (unstranded), 1 (stranded) and 2 (reversely stranded). 0 by default.

  - **t** (str, required) -- Specify the feature type. Only rows which have the matched feature type in the provided GTF annotation file will be included for read counting. 'exon' by default.


**Required tools:** feature_counts

**CPU Cores:** 4

.. index:: filter_gtf

filter_gtf
==========



    custom script to filter merged genocde gtf from cufflinks or stringtie by classcode

**Input Connection**
  - **in/assembling**

**Output Connection**
  - **out/assembling**
  - **out/log_stderr**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      filter_gtf [style=filled, fillcolor="#fce94f"];
      in_0 [label="assembling"];
      in_0 -> filter_gtf;
      out_1 [label="assembling"];
      filter_gtf -> out_1;
      out_2 [label="log_stderr"];
      filter_gtf -> out_2;
   }

**Options:**
  - **class-code-only-in-transcript-feature** (bool, required) -- Removes transcripts without strand specifity

  - **class-list-keep** (str, optional) -- class codes to be kept possible '=,c,j,e,i,o,p,r,u,x,s,.'

  - **keep-by-class** (bool, required) -- "keep gtf if any class is found in class_code field, requieres class-list-keep

  - **remove-by-field-match** (str, optional) -- select gft field like gene_id, gene_name which will match against --string

  - **remove-unstranded** (bool, required) -- Removes transcripts without strand specifity

  - **remove-unwanted-chr** (bool, required) -- keeps chr1 ..2 and chrX, chrY, chrMT

  - **string** (str, optional) -- string to match in gtf field gene_name for discarding


**Required tools:** filter_gtf

**CPU Cores:** 1

.. index:: fix_cutadapt

fix_cutadapt
============



    This step takes FASTQ data and removes both reads of a paired-end read, if
    one of them has been completely removed by cutadapt (or any other software).

**Input Connection**
  - **in/second_read**
  - **in/first_read**

**Output Connection**
  - **out/second_read**
  - **out/first_read**


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
  - **dd-blocksize** (str, optional)    - default value: 2M

  - **pigz-blocksize** (str, optional)    - default value: 2048


**Required tools:** cat (coreutils), dd (coreutils), mkfifo (coreutils), pigz

**CPU Cores:** 4

.. index:: fusioncatcher

fusioncatcher
=============



    FusionCatcher is a tool to discover gene fusions
    in human paired-end RNA-Seq data.

    Paper:
    https://www.biorxiv.org/content/early/2014/11/19/011650

    Manual including required folder structure and typical usage:
    https://github.com/ndaniel/fusioncatcher/blob/master/doc/manual.md


**Input Connection**
  - **in/second_read**
  - **in/first_read**

**Output Connection**
  - **out/tar_archive**
  - **out/log_stdout**
  - **out/log_stderr**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      fusioncatcher [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> fusioncatcher;
      in_1 [label="second_read"];
      in_1 -> fusioncatcher;
      out_2 [label="log_stderr"];
      fusioncatcher -> out_2;
      out_3 [label="log_stdout"];
      fusioncatcher -> out_3;
      out_4 [label="tar_archive"];
      fusioncatcher -> out_4;
   }

**Options:**
  - **cores** (str, required)    - default value: 6

  - **extract-buffer-size** (str, optional)
  - **index** (str, required) -- Path to index folder

  - **keep-unmapped-read** (bool, optional)
  - **skip-filter-adapter** (bool, optional)

**Required tools:** fusioncatcher, mkdir (coreutils), rm (coreutils), tar

**CPU Cores:** 6

.. index:: gffcompare_single

gffcompare_single
=================



    gffcompare [-r <reference_mrna.gtf> [-R]] [-G] [-T] [-V] [-s <seq_path>]
        [-o <outprefix>] [-p <cprefix>]
        {-i <input_gtf_list> | <input1.gtf> [<input2.gtf> .. <inputN.gtf>]}

     GffCompare provides classification and reference annotation mapping and
     matching statistics for RNA-Seq assemblies (transfrags) or other generic
     GFF/GTF files.
     GffCompare also clusters and tracks transcripts across multiple GFF/GTF
     files (samples), writing matching transcripts (identical intron chains) into
     <outprefix>.tracking, and a GTF file <outprefix>.combined.gtf which
     contains a nonredundant set of transcripts across all input files (with
     a single representative transfrag chosen for each clique of matching transfrags
     across samples).

**Input Connection**
  - **in/assembling**

**Output Connection**
  - **out/loci**
  - **out/combined**
  - **out/stats**
  - **out/log_stderr**
  - **out/tracking**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      gffcompare_single [style=filled, fillcolor="#fce94f"];
      in_0 [label="assembling"];
      in_0 -> gffcompare_single;
      out_1 [label="combined"];
      gffcompare_single -> out_1;
      out_2 [label="loci"];
      gffcompare_single -> out_2;
      out_3 [label="log_stderr"];
      gffcompare_single -> out_3;
      out_4 [label="stats"];
      gffcompare_single -> out_4;
      out_5 [label="tracking"];
      gffcompare_single -> out_5;
   }

**Options:**
  - **C** (bool, optional) -- include the "contained" transcripts in the .combined.gtf file

  - **F** (bool, optional) -- do not discard intron-redundant transfrags if they share the 5' end (if they differ only at the 3' end)

  - **G** (bool, optional) -- generic GFF input file(s): do not assume Cufflinks/Stringtie GTF input, (do not discard intron-redundant transfrags)

  - **M** (bool, optional) -- discard (ignore) single-exon transfrags and reference transcripts

  - **N** (bool, optional) -- discard (ignore) single-exon reference transcripts

  - **Q** (bool, optional) -- for -r option, consider only the input transcripts that overlap any of the reference transcripts (Precision correction); (Warning: this will discard all "novel" loci!)

  - **R** (bool, optional) -- for -r option, consider only the reference transcripts that overlap any of the input transfrags (Sn correction)

  - **T** (bool, optional) -- do not generate .tmap and .refmap files for each input file

  - **d** (int, optional) -- max. distance (range) for grouping transcript start sites (100)

  - **e** (int, optional) -- max. distance (range) allowed from free ends of terminal exons of reference transcripts when assessing exon accuracy (100)

  - **i** (str, optional) -- rrovide a text file with a list of (query) GTF files to process instead of expecting them as command line arguments useful when a large number of GTF files should be processed)

  - **p** (str, optional) -- rthe name prefix to use for consensus transcripts in the <outprefix>.combined.gtf file (default: 'TCONS')

  - **r** (str, optional) -- reference annotation file (GTF/GFF)

  - **s** (str, optional) -- path to genome sequences (optional); this can be either a multi-FASTA file or a directory containing single-fasta files (one for each contig); repeats must be soft-masked (lower case) in order to be able to classify transfrags as repeats


**Required tools:** gffcompare, mkdir (coreutils), mv (coreutils)

**CPU Cores:** 2

.. index:: gffread_extract_transcripts

gffread_extract_transcripts
===========================



    extract transcripts from gtf
    http://ccb.jhu.edu/software/stringtie/gff.shtml
    write a fasta file with spliced exons for each GFF transcript
    gffread -w transcripts.fa -g /path/to/genome.fa transcripts.gtf

**Input Connection**
  - **in/fasta**
  - **in/anno**

**Output Connection**
  - **out/fasta**
  - **out/log_stdout**
  - **out/log_stderr**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      gffread_extract_transcripts [style=filled, fillcolor="#fce94f"];
      in_0 [label="anno"];
      in_0 -> gffread_extract_transcripts;
      in_1 [label="fasta"];
      in_1 -> gffread_extract_transcripts;
      out_2 [label="fasta"];
      gffread_extract_transcripts -> out_2;
      out_3 [label="log_stderr"];
      gffread_extract_transcripts -> out_3;
      out_4 [label="log_stdout"];
      gffread_extract_transcripts -> out_4;
   }

**Options:**
  - **gtf** (str, optional) -- path to gtf file

  - **output-fasta-name** (str, required) -- name of the outputfile trancriptom myfasta.fa


**Required tools:** gffread

**CPU Cores:** 1

.. index:: gsnap

gsnap
=====




**Input Connection**
  - **in/second_read**
  - **in/first_read**

**Output Connection**
  - **out/alignments**
  - **out/log_stdout**
  - **out/log_stderr**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      gsnap [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> gsnap;
      in_1 [label="second_read"];
      in_1 -> gsnap;
      out_2 [label="alignments"];
      gsnap -> out_2;
      out_3 [label="log_stderr"];
      gsnap -> out_3;
      out_4 [label="log_stdout"];
      gsnap -> out_4;
   }

**Options:**
  - **D** (str, required) -- Genome directory

  - **cores** (int, optional) -- workaround to specify cores for grid engine and threads ie
    - default value: 1

  - **d** (str, required) -- Genome database

  - **t** (int, optional) -- Number of worker threads
    - default value: 1


**Required tools:** gsnap

**CPU Cores:** 1

.. index:: hisat2

hisat2
======



    HISAT2 is a fast and sensitive alignment program for mapping
    next-generation sequencing reads (both DNA and RNA) to a population of
    human genomes (as well as to a single reference genome).

    https://ccb.jhu.edu/software/hisat2/index.shtml
    must be version 2.1 or higher
    metrics and summary file are automatically produced

**Input Connection**
  - **in/second_read** (optional)
  - **in/first_read**

**Output Connection**
  - **out/metrics**
  - **out/unaligned** (optional) Format: **fastq.gz** - Unpaired reads that didn't align.
  - **out/summary**
  - **out/alignments**
  - **out/aligned** (optional) Format: **fastq.gz** - Unpaired reads that aligned.
  - **out/log_stderr**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      hisat2 [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> hisat2;
      in_1 [label="second_read", style=filled, fillcolor="#a7a7a7"];
      in_1 -> hisat2;
      out_2 [label="aligned", style=filled, fillcolor="#a7a7a7"];
      hisat2 -> out_2;
      out_3 [label="alignments"];
      hisat2 -> out_3;
      out_4 [label="log_stderr"];
      hisat2 -> out_4;
      out_5 [label="metrics"];
      hisat2 -> out_5;
      out_6 [label="summary"];
      hisat2 -> out_6;
      out_7 [label="unaligned", style=filled, fillcolor="#a7a7a7"];
      hisat2 -> out_7;
   }

**Options:**
  - **add-chrname** (bool, optional) -- Add 'chr' to reference names in alignment (e.g., 18 to chr18)

  - **al-gz** (bool, optional) -- write unpaired reads that aligned to gzip compress output connection "out/aligned"

  - **c** (bool, optional) -- <m1>, <m2>, <r> are sequences themselves, not files

  - **cores** (int, required)    - default value: 12

  - **dta** (bool, optional) -- Reports alignments tailored for transcript assemblers

  - **f** (bool, optional) -- query input files are (multi-)FASTA .fa/.mfa

  - **ff** (bool, optional) -- -1, -2 mates align fw/rev, rev/fw, fw/fw (--fr)

  - **fr** (bool, optional) -- -1, -2 mates align fw/rev, rev/fw, fw/fw (--fr)

  - **ignore-quals** (bool, optional) -- treat all quality values as 30 on Phred scale (off)

  - **index** (str, required) -- Path to hisat2 index (not containing file suffixes).

  - **int-quals** (bool, optional) -- qualities encoded as space-delimited integers

  - **k** (int, optional) -- report up to <int> alns per read; MAPQ not meaningful

  - **known-splicesite-infile** (str, optional) -- provide a list of known splice sites

  - **library_type** (str, optional) -- -1, -2 mates align fr (fw/rev), rf (rev/fw), ff (fw/fw) (default fr).
    - possible values: 'fr', 'rf', 'ff'


  - **ma** (str, optional) -- match bonus (0 for --end-to-end, 2 for --local)

  - **max-intronlen** (str, optional) -- maximum intron length (500000)

  - **maxins** (int, optional) -- maximum fragment length (500), only valid with --no-spliced-alignment

  - **min-intronlen** (str, optional) -- minimum intron length (20)

  - **minins** (int, optional) -- minimum fragment length (0), only valid with --no-spliced-alignment

  - **mm** (bool, optional) -- use memory-mapped I/O for index; many 'hisat2's can share

  - **mp** (str, optional) -- max and min penalties for mismatch; lower qual = lower penalty <2,6>

  - **n-ceil** (str, optional) -- func for max # non-A/C/G/Ts permitted in aln (L,0,0.15)

  - **new-summary** (bool, optional) -- print alignment summary in a new style, which is more machine-friendly

  - **no-discordant** (bool, optional) -- suppress discordant alignments for paired reads

  - **no-head** (bool, optional) -- supppress header lines, i.e. lines starting with @

  - **no-mixed** (bool, optional) -- suppress unpaired alignments for paired reads

  - **no-softclip** (bool, optional) -- no soft-clipping

  - **no-spliced-alignment** (bool, optional) -- disable spliced alignment

  - **no-sq** (bool, optional) -- supppress @SQ header lines

  - **no-temp-splicesite** (bool, optional) -- disable the use of splice sites found

  - **nofw** (bool, optional) -- do not align forward (original) version of read (off)

  - **non-deterministic** (bool, optional) -- seed rand. gen. arbitrarily instead of using read attributes

  - **norc** (bool, optional) -- do not align reverse-complement version of read (off)

  - **novel-splicesite-infile** (str, optional) -- provide a list of novel splice sites

  - **novel-splicesite-outfile** (str, optional) -- report a list of splice sites

  - **np** (str, optional) -- penalty for non-A/C/G/Ts in read/ref (1)

  - **offrate** (int, optional) -- override offrate of index; must be >= index's offrate

  - **omit-sec-seq** (bool, optional) -- put '*' in SEQ and QUAL fields for secondary alignments

  - **pen-canintronlen** (str, optional) -- penalty for long introns (G,-8,1) with canonical splice sites

  - **pen-cansplice** (str, optional) -- penalty for a canonical splice site (0)

  - **pen-noncanintronlen** (str, optional) -- penalty for long introns (G,-8,1) with noncanonical splice sites

  - **pen-noncansplice** (str, optional) -- penalty for a non-canonical splice site (12)

  - **phred33** (bool, optional) -- qualities are Phred+33 (default)

  - **phred64** (bool, optional) -- qualities are Phred+64

  - **q** (bool, optional) -- query input files are FASTQ .fq/.fastq (default)

  - **qc-filter** (bool, optional) -- filter out reads that are bad according to QSEQ filter

  - **qseq** (bool, optional) -- query input files are in Illumina's qseq format

  - **quiet** (bool, optional) -- print nothing to stderr except serious errors

  - **r** (bool, optional) -- query input files are raw one-sequence-per-line

  - **rdg** (str, optional) -- read gap open, extend penalties (5,3)

  - **remove-chrname** (bool, optional) -- Remove 'chr' from reference names in alignment (e.g., chr18 to 18)

  - **reorder** (bool, optional) -- force SAM output order to match order of input reads

  - **rf** (bool, optional) -- -1, -2 mates align fw/rev, rev/fw, fw/fw (--fr)

  - **rfg** (str, optional) -- reference gap open, extend penalties (5,3)

  - **rg** (str, optional) -- add <text> ('lab:value') to @RG line of SAM header. (Note: @RG line only printed when --rg-id is set.)

  - **rg-id** (str, optional) -- puts sample name in rg

  - **rna-strandness** (str, required) -- Specify strand-specific information (unstranded); paired and are extended F->FR, R->RF
    - possible values: 'R', 'F', 'U'


  - **score-min** (str, optional) -- min acceptable alignment score w/r/t read length (G,20,8 for local, L,-0.6,-0.6 for end-to-end)

  - **seed** (int, optional) -- seed for random number generator (0)

  - **skip** (int, optional) -- skip the first <int> reads/pairs in the input (none)

  - **sp** (str, optional) -- max and min penalties for soft-clipping; lower qual = lower penalty <1,2>

  - **tmo** (bool, optional) -- Reports only those alignments within known transcriptome

  - **trim3** (int, optional) -- trim <int> bases from 3'/right end of reads (0)

  - **trim5** (int, optional) -- trim <int> bases from 5'/left end of reads (0)

  - **un-gz** (bool, optional) -- write unpaired reads that didn't align to gzip compress output connection "out/unaligned"

  - **upto** (int, optional) -- stop after first <int> reads/pairs (no limit)


**Required tools:** hisat2, pigz

**CPU Cores:** 12

.. index:: htseq_count

htseq_count
===========



    The htseq-count script counts the number of reads overlapping a feature.
    Input needs to be a file with aligned sequencing reads and a list of genomic
    features. For more information see:

    http://www-huber.embl.de/users/anders/HTSeq/doc/count.html

**Input Connection**
  - **in/features** (optional) Format: **gtf** - reference assembly
  - **in/alignments**

**Output Connection**
  - **out/counts**


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
      in_1 [label="features", style=filled, fillcolor="#a7a7a7"];
      in_1 -> htseq_count;
      out_2 [label="counts"];
      htseq_count -> out_2;
   }

**Options:**
  - **a** (int, optional)
  - **dd-blocksize** (str, optional)    - default value: 2M

  - **feature-file** (str, optional)
  - **idattr** (str, optional)    - default value: gene_id

  - **mode** (str, optional)    - default value: union
    - possible values: 'union', 'intersection-strict', 'intersection-nonempty'


  - **order** (str, required)    - possible values: 'name', 'pos'


  - **pigz-blocksize** (str, optional)    - default value: 2048

  - **stranded** (str, required)    - possible values: 'yes', 'no', 'reverse'


  - **threads** (int, optional) -- start <n> threads (default:2)
    - default value: 2

  - **type** (str, optional)    - default value: exon


**Required tools:** dd (coreutils), htseq-count, pigz, samtools

**CPU Cores:** 2

.. index:: identify_adapters

identify_adapters
=================



    Uses AdapterRemoval to identify adapter sequences from paired read data.

    AdapterRemoval (ver. 2.1.7)
    This program searches for and removes remnant adapter sequences from
    your read data.  The program can analyze both single end and paired end
    data.  For detailed explanation of the parameters, please refer to the
    man page.  For comments, suggestions  and feedback please contact Stinus
    Lindgreen (stinus@binf.ku.dk) and Mikkel Schubert (MikkelSch@gmail.com).
    If you use the program, please cite the paper:
    Schubert, Lindgreen, and Orlando (2016). AdapterRemoval v2: rapid
    adapter trimming, identification, and read merging.
    BMC Research Notes, 12;9(1):88.
    http://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-016-1900-2
    "Pipeline specific "input and output expected to be gzipped"

**Input Connection**
  - **in/second_read**
  - **in/first_read**

**Output Connection**
  - **out/log_stdout**
  - **out/log_stderr**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      identify_adapters [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> identify_adapters;
      in_1 [label="second_read"];
      in_1 -> identify_adapters;
      out_2 [label="log_stderr"];
      identify_adapters -> out_2;
      out_3 [label="log_stdout"];
      identify_adapters -> out_3;
   }

**Required tools:** adapterremoval

**CPU Cores:** 1

.. index:: kallisto

kallisto
========



    Kallisto



**Input Connection**
  - **in/second_read** (optional)
  - **in/kallisto-index** (optional)
  - **in/first_read**

**Output Connection**
  - **out/run_info.json**
  - **out/abundance.h5**
  - **out/log_stdout**
  - **out/log_stderr**
  - **out/abundance.tsv**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      kallisto [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> kallisto;
      in_1 [label="kallisto-index", style=filled, fillcolor="#a7a7a7"];
      in_1 -> kallisto;
      in_2 [label="second_read", style=filled, fillcolor="#a7a7a7"];
      in_2 -> kallisto;
      out_3 [label="abundance.h5"];
      kallisto -> out_3;
      out_4 [label="abundance.tsv"];
      kallisto -> out_4;
      out_5 [label="log_stderr"];
      kallisto -> out_5;
      out_6 [label="log_stdout"];
      kallisto -> out_6;
      out_7 [label="run_info.json"];
      kallisto -> out_7;
   }

**Options:**
  - **bias** (bool, optional) -- Perform sequence based bias correction

  - **bootstrap-samples** (int, optional) -- Number of bootstrap samples (default: 0)

  - **cores** (int, optional) -- workaround to specify cores for grid engine and threads ie
    - default value: 1

  - **fr-stranded** (bool, optional) -- Strand specific reads, first read forward

  - **fragment-length** (int, optional) -- Estimated average fragment length

  - **index** (str, optional) -- Filename for the kallisto index to be used for quantification

  - **rf-stranded** (bool, optional) -- Strand specific reads, first read reverse

  - **sd** (int, optional) -- Estimated standard deviation of fragment length

  - **seed** (int, optional) -- Seed for the bootstrap sampling

  - **single** (bool, optional) -- Quantify single-end reads

  - **single-overhang** (bool, optional) -- Include reads where unobserved rest of fragment is predicted to lie outside a transcript


**Required tools:** kallisto

**CPU Cores:** 1

.. index:: kallisto_fusion

kallisto_fusion
===============



    Kallisto with fusion detection option



**Input Connection**
  - **in/second_read**
  - **in/kallisto-index**
  - **in/first_read**

**Output Connection**
  - **out/run_info.json**
  - **out/abundance.h5**
  - **out/fusion.txt**
  - **out/log_stdout**
  - **out/log_stderr**
  - **out/abundance.tsv**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      kallisto_fusion [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> kallisto_fusion;
      in_1 [label="kallisto-index"];
      in_1 -> kallisto_fusion;
      in_2 [label="second_read"];
      in_2 -> kallisto_fusion;
      out_3 [label="abundance.h5"];
      kallisto_fusion -> out_3;
      out_4 [label="abundance.tsv"];
      kallisto_fusion -> out_4;
      out_5 [label="fusion.txt"];
      kallisto_fusion -> out_5;
      out_6 [label="log_stderr"];
      kallisto_fusion -> out_6;
      out_7 [label="log_stdout"];
      kallisto_fusion -> out_7;
      out_8 [label="run_info.json"];
      kallisto_fusion -> out_8;
   }

**Options:**
  - **bias** (bool, optional) -- Perform sequence based bias correction

  - **bootstrap-samples** (int, optional) -- Number of bootstrap samples (default: 0)

  - **cores** (int, optional) -- workaround to specify cores for grid engine and threads ie
    - default value: 1

  - **fr-stranded** (bool, optional) -- Strand specific reads, first read forward

  - **fragment-length** (int, optional) -- Estimated average fragment length

  - **index** (str, required) -- Filename for the kallisto index to be used for quantification

  - **rf-stranded** (bool, optional) -- Strand specific reads, first read reverse

  - **sd** (int, optional) -- Estimated standard deviation of fragment length

  - **seed** (int, optional) -- Seed for the bootstrap sampling

  - **single** (int, optional) -- Quantify single-end reads

  - **single-overhang** (int, optional) -- Include reads where unobserved rest of fragment is predicted to lie outside a transcript


**Required tools:** kallisto

**CPU Cores:** 1

.. index:: kallisto_index

kallisto_index
==============



    https://pachterlab.github.io/kallisto/manual
    Builds a kallisto index
    Usage: kallisto index [arguments] FASTA-files

**Input Connection**
  - **in/fasta**

**Output Connection**
  - **out/log_stdout**
  - **out/log_stderr**
  - **out/kallisto-index**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      kallisto_index [style=filled, fillcolor="#fce94f"];
      in_0 [label="fasta"];
      in_0 -> kallisto_index;
      out_1 [label="kallisto-index"];
      kallisto_index -> out_1;
      out_2 [label="log_stderr"];
      kallisto_index -> out_2;
      out_3 [label="log_stdout"];
      kallisto_index -> out_3;
   }

**Options:**
  - **index** (str, required) -- Filename for the kallisto index

  - **kmer-size** (int, optional) -- k-mer (odd) length (default: 31, max value: 31)
    - default value: 31

  - **make-unique** (bool, optional) -- Replace repeated target names with unique names


**Required tools:** kallisto

**CPU Cores:** 1

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

    typical command line for single-end data::

        macs2 callpeak --treatment <aligned-reads> [--control <aligned-reads>]
                       --name <run-id> --gsize 2.7e9

**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/model**
  - **out/broadpeaks-xls**
  - **out/gappedpeaks**
  - **out/summits**
  - **out/narrowpeaks-xls**
  - **out/narrowpeaks**
  - **out/diagnosis**
  - **out/log**
  - **out/broadpeaks**


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
  - **bdg** (bool, optional) -- If this flag is on, MACS will store the fragment pileup, control lambda, -log10pvalue and -log10qvalue scores in bedGraph files. The bedGraph files will be stored in current directory named NAME+'_treat_pileup.bdg' for treatment data, NAME+'_control_lambda.bdg' for local lambda values from control, NAME+'_treat_pvalue.bdg' for Poisson pvalue scores (in -log10(pvalue) form), and NAME+'_treat_qvalue.bdg' for q-value scores from Benjamini-Hochberg-Yekutieli procedure <http://en.wikipedia.org/wiki/False_discovery_rate#Dependent_tests>

  - **broad** (bool, optional) -- When this flag is on, MACS will try to composite broad regions in BED12 ( a gene-model-like format ) by putting nearby highly enriched regions into a broad region with loose cutoff. The broad region is controlled by another cutoff through --broad-cutoff. The maximum length of broad region length is 4 times of d from MACS. DEFAULT: False

  - **broad-cutoff** (float, optional) -- Cutoff for broad region. This option is not available unless --broad is set. If -p is set, this is a pvalue cutoff, otherwise, it's a qvalue cutoff. DEFAULT: 0.1

  - **buffer-size** (int, optional) -- LEGACY option.

  - **bw** (int, optional) -- The band width which is used to scan the genome ONLY for model building. You can set this parameter as the sonication fragment size expected from wet experiment. The previous side effect on the peak detection process has been removed. So this parameter only affects the model building.

  - **call-summits** (bool, optional) -- MACS will now reanalyze the shape of signal profile (p or q-score depending on cutoff setting) to deconvolve subpeaks within each peak called from general procedure. It's highly recommended to detect adjacent binding events. While used, the output subpeaks of a big peak region will have the same peak boundaries, and different scores and peak summit positions.

  - **control** (dict, required) -- Defines the controls and correspondent treatments in a YAML hash. Hash keys are the run IDs of the control datasets and hash values are the run IDs of the treatment datasets.

  - **down-sample** (bool, optional) -- When set, random sampling method will scale down the bigger sample. By default, MACS uses linear scaling. This option will make the results unstable and irreproducible since each time, random reads would be selected, especially the numbers (pileup, pvalue, qvalue) would change. Consider to use 'randsample' script before MACS2 runs instead.

  - **extsize** (int, optional) -- While '--nomodel' is set, MACS uses this parameter to extend reads in 5'->3' direction to fix-sized fragments. For example, if the size of binding region for your transcription factor is 200 bp, and you want to bypass the model building by MACS, this parameter can be set as 200. This option is only valid when --nomodel is set or when MACS fails to build model and --fix-bimodal is on.

  - **fix-bimodal** (bool, optional) -- Whether turn on the auto paired-peak model process. If it's set, when MACS failed to build paired model, it will use the nomodel settings, the '--extsize' parameter to extend each tags. If set, MACS will be terminated if paried-peak model is failed.

  - **format** (str, required) -- Format of tag file, can be 'ELAND', 'BED', 'ELANDMULTI', 'ELANDEXPORT', 'ELANDMULTIPET' (for pair-end tags), 'SAM', 'BAM', 'BOWTIE', 'BAMPE' or 'BEDPE'. Default is 'AUTO' which will allow MACS to decide the format automatically. 'AUTO' is also useful when you combine different formats of files. Note that MACS can't detect 'BAMPE' or 'BEDPE' format with 'AUTO', and you have to implicitly specify the format for 'BAMPE' and 'BEDPE'. For more information about the formats see https://github.com/taoliu/MACS/
    - default value: AUTO
    - possible values: 'AUTO', 'ELAND', 'ELANDMULTI', 'ELANDMULTIPET', 'ELANDEXPORT', 'BED', 'BEDPE', 'SAM', 'BAM', 'BAMPE', 'BOWTIE'


  - **gsize** (str, required) -- PLEASE assign this parameter to fit your needs! It's the mappable genome size or effective genome size which is defined as the genome size which can be sequenced. Because of the repetitive features on the chromsomes, the actual mappable genome size will be smaller than the original size, about 90% or 70% of the genome size. The default hs -- 2.7e9 is recommended for UCSC human hg18 assembly. Here are all precompiled parameters for effective genome size: hs:2.7e9; mm:1.87e9; ce:9e7; dm:1.2e8
    - default value: 2.7e9

  - **keep-dup** (int, optional) -- It controls the MACS behavior towards duplicate tags at the exact same location -- the same coordination and the same strand. The default 'auto' option makes MACS calculate the maximum tags at the exact same location based on binomal distribution using 1e-5 as pvalue cutoff; and the 'all' option keeps every tags. If an integer is given, at most this number of tags will be kept at the same location. The default is to keep one tag at the same location. Default: 1

  - **llocal** (str, optional) -- 'slocal' and 'llocal' control which two levels of regions will be checked around the peak regions to calculate the maximum lambda as local lambda. By default, MACS considers 1000bp for small local region(--slocal), and 10000bps for large local region(--llocal) which captures the bias from a long range effect like an open chromatin domain. You can tweak these according to your project. Remember that if the region is set too small, a sharp spike in the input data may kill the significant peak.

  - **mfold** (str, optional) -- This parameter is used to select the regions within MFOLD range of high-confidence enrichment ratio against background to build model. The regions must be lower than upper limit, and higher than the lower limit of fold enrichment. DEFAULT:5,50 means using all regions not too low (>5) and not too high (<50) to build paired-peaks model. If MACS can not find more than 100 regions to build model, it will use the --extsize parameter to continue the peak detection ONLY if --fix-bimodal is set.

  - **nolambda** (bool, optional) -- With this flag on, MACS will use the background lambda as local lambda. This means MACS will not consider the local bias at peak candidate regions.

  - **nomodel** (bool, optional) -- While on, MACS will bypass building the shifting model.

  - **pvalue** (float, optional) -- The pvalue cutoff. If 'pvalue' is specified, MACS2 will use pvalue instead of qvalue.

  - **qvalue** (float, optional) -- The qvalue (minimum FDR) cutoff to call significant regions. Default is 0.05. For broad marks, you can try 0.05 as cutoff. Q-values are calculated from p-values using Benjamini-Hochberg procedure.

  - **read-length** (int, optional) -- LEGACY option.

  - **shift** (int, optional)
  - **slocal** (str, optional) -- 'slocal' and 'llocal' control which two levels of regions will be checked around the peak regions to calculate the maximum lambda as local lambda. By default, MACS considers 1000bp for small local region(--slocal), and 10000bps for large local region(--llocal) which captures the bias from a long range effect like an open chromatin domain. You can tweak these according to your project. Remember that if the region is set too small, a sharp spike in the input data may kill the significant peak.

  - **to-large** (bool, optional) -- When set, linearly scale the smaller dataset to the same depth as larger dataset, by default, the larger dataset will be scaled towards the smaller dataset. Beware, to scale up small data would cause more false positives.

  - **tsize** (int, optional) -- The size of sequencing tags. If you don't specify it, MACS will try to use the first 10 sequences from your input treatment file to determine the tag size. Specifying it will override the automatically determined tag size.

  - **verbose** (int, optional) -- If you don't want to see any message during the running of MACS, set it to 0. But the CRITICAL messages will never be hidden. If you want to see rich information like how many peaks are called for every chromosome, you can set it to 3 or larger than 3.
    - possible values: '0', '1', '2', '3'



**Required tools:** macs2, mkdir (coreutils), mv (coreutils)

**CPU Cores:** 4

.. index:: merge_assembly

merge_assembly
==============



    This step merges a single gtf file that has been produced by a previous step
    (e.g. cuffmerge or cuffcompare) with a reference annotation. No lines
    are discarded. The two files are simply concatenated and subsequently
    sorted by position.


**Input Connection**
  - **in/features**

**Output Connection**
  - **out/features**
  - **out/log_stderr**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      merge_assembly [style=filled, fillcolor="#fce94f"];
      in_0 [label="features"];
      in_0 -> merge_assembly;
      out_1 [label="features"];
      merge_assembly -> out_1;
      out_2 [label="log_stderr"];
      merge_assembly -> out_2;
   }

**Options:**
  - **reference** (str, required) -- The reference annotation file that should be merged.

  - **temp-sort-dir** (str, required) -- Intermediate sort files are stored into this directory. Note that this directory needs to be present before running this step.


**Required tools:** cat (coreutils), sort (coreutils)

**CPU Cores:** 1

.. index:: merge_fasta_files

merge_fasta_files
=================



    This step concatenates all .fasta(.gz) files that belong to a certain
    sample.

**Input Connection**
  - **in/sequence**

**Output Connection**
  - **out/sequence**


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
  - **compress-output** (bool, optional) -- Produce gzipped output.
    - default value: True

  - **dd-blocksize** (str, optional)    - default value: 2M

  - **merge-all-runs** (bool, optional) -- Merge sequences from all runs.

  - **output-fasta-basename** (str, optional) -- Name used as prefix for FASTA output.

  - **pigz-blocksize** (str, optional)    - default value: 2048


**Required tools:** cat (coreutils), dd (coreutils), mkfifo (coreutils), pigz

**CPU Cores:** 4

.. index:: merge_fastq_files

merge_fastq_files
=================



    This step concatenates all .fastq(.gz) files belonging to a certain sample.
    First and second read files are merged separately. The output files are
    gzipped.

**Input Connection**
  - **in/second_read**
  - **in/first_read**

**Output Connection**
  - **out/second_read**
  - **out/first_read**


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
  - **dd-blocksize** (str, optional)    - default value: 2M

  - **pigz-blocksize** (str, optional)    - default value: 2048


**Required tools:** cat (coreutils), dd (coreutils), mkfifo (coreutils), pigz

**CPU Cores:** 4

.. index:: merge_fastx_files

merge_fastx_files
=================



    This step merges all .fastq/a(.gz) files belonging to a certain sample.
    First and second read files are merged separately. The output files are
    gzipped.

**Input Connection**
  - **in/second_read**
  - **in/first_read**

**Output Connection**
  - **out/first_read**
  - **out/second_read**
  - **out/report**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      merge_fastx_files [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> merge_fastx_files;
      in_1 [label="second_read"];
      in_1 -> merge_fastx_files;
      out_2 [label="first_read"];
      merge_fastx_files -> out_2;
      out_3 [label="report"];
      merge_fastx_files -> out_3;
      out_4 [label="second_read"];
      merge_fastx_files -> out_4;
   }

**Required tools:** echo, pigz

**CPU Cores:** 4

.. index:: merge_genecounts

merge_genecounts
================





**Input Connection**
  - **in/counts**

**Output Connection**
  - **out/merged_counts**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      merge_genecounts [style=filled, fillcolor="#fce94f"];
      in_0 [label="counts"];
      in_0 -> merge_genecounts;
      out_1 [label="merged_counts"];
      merge_genecounts -> out_1;
   }

**Options:**
  - **cores** (int, optional) -- workaround to specify cores for grid engine and threads ie
    - default value: 1

  - **t** (str, required) -- tool name (htseq_count: htc, featureCounts: fc)


**CPU Cores:** 1

.. index:: merge_numpy_zip_arrays

merge_numpy_zip_arrays
======================



    This step can be used to concatenate multiple zipped Numpy arrays which are
    the output of deepTools multiBamSummary subcommand.

    Usage example::

        merge_numpy_arrays.py [-h] [file-1.npz file-2.npz ... file-n.npz]  files


**Input Connection**
  - **in/read-coverage**

**Output Connection**
  - **out/read-coverage**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      merge_numpy_zip_arrays [style=filled, fillcolor="#fce94f"];
      in_0 [label="read-coverage"];
      in_0 -> merge_numpy_zip_arrays;
      out_1 [label="read-coverage"];
      merge_numpy_zip_arrays -> out_1;
   }

**Required tools:** merge_numpy_arrays.py

**CPU Cores:** 10

.. index:: pear

pear
====


**Input Connection**
  - **in/second_read**
  - **in/first_read**

**Output Connection**
  - **out/unassembled.forward**
  - **out/unassembled.reverse**
  - **out/discarded**
  - **out/assembled**
  - **out/log**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      pear [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> pear;
      in_1 [label="second_read"];
      in_1 -> pear;
      out_2 [label="assembled"];
      pear -> out_2;
      out_3 [label="discarded"];
      pear -> out_3;
      out_4 [label="log"];
      pear -> out_4;
      out_5 [label="unassembled.forward"];
      pear -> out_5;
      out_6 [label="unassembled.reverse"];
      pear -> out_6;
   }

**Options:**
  - **cap** (int, optional) -- Specify the upper bound for the resulting quality score.
    - default value: 40
    - possible values: '0', '33', '34', '35', '36', '37', '38', '39', '40', '41', '42', '43', '44', '45', '46', '47', '48', '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', '61', '62', '63', '64', '65', '66', '67', '68', '69', '70', '71', '72', '73', '74', '75', '76', '77', '78', '79', '80', '81', '82', '83', '84', '85', '86', '87', '88', '89', '90', '91', '92', '93', '94', '95', '96', '97', '98', '99', '100', '101', '102', '103', '104', '105', '106', '107', '108', '109', '110', '111', '112', '113', '114', '115', '116', '117', '118', '119', '120', '121', '122', '123', '124', '125', '126', '127', '128', '129', '130', '131', '132', '133', '134', '135', '136', '137', '138', '139', '140', '141', '142', '143', '144', '145', '146', '147', '148', '149', '150', '151', '152', '153', '154', '155', '156', '157'


  - **empirical-freqs** (bool, optional) -- Disable empirical basefrequencies.

  - **m** (int, optional) -- Specify the maximum possible length of the assembled sequences.
    - default value: 300

  - **max-uncalled-base** (float, optional) -- Specify the maximal proportion of uncalled bases in a read.
    - default value: 1.0
    - possible values: '[0.0..1.0]'


  - **memory** (str, optional) -- Specify the amount of memory to be used.

  - **min-overlap** (int, optional) -- Specify the minimum overlap size.
    - default value: 10

  - **min-trim-length** (int, optional) -- Specify the minimum length of reads after trimming the low quality part.
    - default value: 1

  - **n** (int, optional) -- Specify the minimum possible length of the assembled sequences.
    - default value: 50

  - **nbase** (bool, optional) -- When resolving a mismatching base-pair out of which non is degenerate, set the merged base to N and use the highest quality score of the two bases.

  - **p-value** (float, optional) -- Specify a p-value for the statistical test.
    - default value: 0.01
    - possible values: '0.0001', '0.001', '0.01', '0.05', '1.0'


  - **phred-base** (int, optional) -- Base PHRED quality score.
    - default value: 33
    - possible values: '33', '64'


  - **q** (int, optional) -- Specify the quality score threshold for trimming the low quality part of a read.
    - possible values: '0', '33', '34', '35', '36', '37', '38', '39', '40', '41', '42', '43', '44', '45', '46', '47', '48', '49', '50', '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', '61', '62', '63', '64', '65', '66', '67', '68', '69', '70', '71', '72', '73', '74', '75', '76', '77', '78', '79', '80', '81', '82', '83', '84', '85', '86', '87', '88', '89', '90', '91', '92', '93', '94', '95', '96', '97', '98', '99', '100', '101', '102', '103', '104', '105', '106', '107', '108', '109', '110', '111', '112', '113', '114', '115', '116', '117', '118', '119', '120', '121', '122', '123', '124', '125', '126', '127', '128', '129', '130', '131', '132', '133', '134', '135', '136', '137', '138', '139', '140', '141', '142', '143', '144', '145', '146', '147', '148', '149', '150', '151', '152', '153', '154', '155', '156', '157'


  - **score-method** (int, optional) -- Specify the scoring method.
    - default value: 2
    - possible values: '1', '2', '3'


  - **test-method** (int, optional) -- Specify the type of statistical tesit.
    - default value: 1
    - possible values: '1', '2'


  - **threads** (int, optional) -- Number of threads to use.


**Required tools:** pear

**CPU Cores:** 1

.. index:: pepr

pepr
====





**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/peaks**
  - **out/parameter**
  - **out/differential_peaks**
  - **out/log**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      pepr [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> pepr;
      out_1 [label="differential_peaks"];
      pepr -> out_1;
      out_2 [label="log"];
      pepr -> out_2;
      out_3 [label="parameter"];
      pepr -> out_3;
      out_4 [label="peaks"];
      pepr -> out_4;
   }

**Options:**
  - **chip_vs_input** (dict, required) -- A YAML dictionary that contains: runID: rep1: [<List of runIDs>] input1: [<List of runIDs>] [ rep2: [<List of runIDs>] input2: [<List of runIDs>] ]rep2 and input2 are optional and will only be used for differential peak calling

  - **diff** (bool, required) -- Tell PePr to perform differential binding analysis or not.

  - **file-format** (str, required) -- Read file format. Currently support bed, sam, bam, sampe (sam paired-end), bampe (bam paired-end)

  - **normalization** (str, optional) -- inter-group, intra-group, scale, or no. Default is intra-group for peak-calling and inter-group for differential binding analysis. PePr is using a modified TMM method to normalize for the difference in IP efficiencies between samples (see the supplementary methods of the paper). It is making an implicit assumption that there is substantial overlap of peaks in every sample. However, it is sometimes not true between groups (for example, between TF ChIP-seq and TF knockout). So for differential binding analysis, switch to intra-group normalization. scale is simply scaling the reads so the total library sizes are the same. no normalization will not do normalization.
    - possible values: 'inter-group', 'intra-group', 'scale', 'no'


  - **peaktype** (str, optional) -- sharp or broad. Default is broad. PePr treats broad peaks (like H3K27me3) and sharp peaks (like most transcriptions factors) slightly different. Specify this option if you know the feature of the peaks.
    - possible values: 'sharp', 'broad'


  - **shiftsize** (str, optional) -- Half the fragment size. The number of bases to shift forward and reverse strand reads toward each other. If not specified by user, PePr will empirically estimate this number from the data for each ChIP sample.

  - **threshold** (float, optional) -- p-value cutoff. Default:1e-5.

  - **windowsize** (int, optional) -- Sliding window size. If not specified by user, PePr will estimate this by calculating the average width of potential peaks. The lower and upper bound for PePr estimate is 100bp and 1000bp. User provided window size is not constrained, but we recommend to stay in this range (100-1000bp).


**Required tools:** mkdir (coreutils), mv (coreutils), pepr, tar

**CPU Cores:** 4

.. index:: pepr_postprocess

pepr_postprocess
================



    Post processing of peaks called by PePr.
    Attention: Filter criteria are hard coded!

**Input Connection**
  - **in/peaks**
  - **in/alignments**

**Output Connection**
  - **out/passed_peaks**
  - **out/peaks**
  - **out/chip**
  - **out/input**
  - **out/failed_peaks**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      pepr_postprocess [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> pepr_postprocess;
      in_1 [label="peaks"];
      in_1 -> pepr_postprocess;
      out_2 [label="chip"];
      pepr_postprocess -> out_2;
      out_3 [label="failed_peaks"];
      pepr_postprocess -> out_3;
      out_4 [label="input"];
      pepr_postprocess -> out_4;
      out_5 [label="passed_peaks"];
      pepr_postprocess -> out_5;
      out_6 [label="peaks"];
      pepr_postprocess -> out_6;
   }

**Options:**
  - **chip_vs_input** (dict, required) -- A YAML dictionary that contains: runID: rep1: [<List of runIDs>] input1: [<List of runIDs>]

  - **file-type** (str, required) -- Read file format. Currently support bed, sam, bam
    - possible values: 'bed', 'sam', 'bam'


  - **narrow-peak-boundary** (bool, optional)
  - **remove-artefacts** (bool, optional)    - default value: True


**Required tools:** ln (coreutils), pepr-postprocess

**CPU Cores:** 4

.. index:: picard_add_replace_read_groups

picard_add_replace_read_groups
==============================



    Replace read groups in a BAM file. This tool enables the user to replace all
    read groups in the INPUT file with a single new read group and assign all
    reads to this read group in the OUTPUT BAM file.

    Documentation::

        https://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups


**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/alignments**


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

.. index:: picard_collectrnaseqmetrics

picard_collectrnaseqmetrics
===========================



    Produces RNA alignment metrics for a SAM or BAM file.

    This tool takes a SAM/BAM file containing the aligned reads from an RNAseq experiment and produces metrics describing the distribution of the bases within the transcripts. It calculates the total numbers and the fractions of nucleotides within specific genomic regions including untranslated regions (UTRs), introns, intergenic sequences (between discrete genes), and peptide-coding sequences (exons). This tool also determines the numbers of bases that pass quality filters that are specific to Illumina data (PF_BASES). For more information please see the corresponding GATK Dictionary entry.

    Other metrics include the median coverage (depth), the ratios of 5 prime /3 prime-biases, and the numbers of reads with the correct/incorrect strand designation. The 5 prime /3 prime-bias results from errors introduced by reverse transcriptase enzymes during library construction, ultimately leading to the over-representation of either the 5 prime or 3 prime ends of transcripts. Please see the CollectRnaSeqMetrics definitions for details on how these biases are calculated.

    The sequence input must be a valid SAM/BAM file containing RNAseq data aligned by an RNAseq-aware genome aligner such a STAR or TopHat. The tool also requires a REF_FLAT file, a tab-delimited file containing information about the location of RNA transcripts, exon start and stop sites, etc. For more information on the REF_FLAT format, see the following description. Build-specific REF_FLAT files can be obtained here.
    Usage example:

    java -jar picard.jar CollectRnaSeqMetrics            I=input.bam            O=output.RNA_Metrics
           REF_FLAT=ref_flat.txt            STRAND=SECOND_READ_TRANSCRIPTION_STRAND            RIBOSOMAL_INTERVALS=ribosomal.interval_list

       https://broadinstitute.github.io/picard/command-line-overview.html#CollectRnaSeqMetrics


**Input Connection**
  - **in/refFlat** (optional)
  - **in/alignments**

**Output Connection**
  - **out/metrics**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      picard_collectrnaseqmetrics [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> picard_collectrnaseqmetrics;
      in_1 [label="refFlat", style=filled, fillcolor="#a7a7a7"];
      in_1 -> picard_collectrnaseqmetrics;
      out_2 [label="metrics"];
      picard_collectrnaseqmetrics -> out_2;
   }

**Options:**
  - **ASSUME_SORTED** (bool, optional)
  - **COMPRESSION_LEVEL** (int, optional) -- Compression level for all compressed files created (e.g. BAM and GELI). Default value: 5. This option can be set to "null" to clear the default value.

  - **CREATE_INDEX** (bool, optional) -- Whether to create a BAM index when writing a coordinate-sorted BAM file. Default value: false. This option can be set to "null" to clear the default value. 

  - **CREATE_MD5_FILE** (bool, optional) -- Whether to create an MD5 digest for any BAM or FASTQ files created. Default value: false. This option can be set to "null" to clear the default value.

  - **GA4GH_CLIENT_SECRETS** (str, optional) -- Google Genomics API client_secrets.json file path. Default value: client_secrets.json. This option can be set to "null" to clear the default value.

  - **IGNORE_SEQUENCE** (str, optional)
  - **MAX_RECORDS_IN_RAM** (int, optional) -- When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed. Default value: 500000. This option can be set to "null" to clear the default value.

  - **METRIC_ACCUMULATION_LEVEL** (str, optional)    - possible values: 'ALL_READS', 'SAMPLE', 'LIBRARY', 'READ_GROUP'


  - **MINIMUM_LENGTH** (int, optional)
  - **QUIET** (bool, optional) -- Whether to suppress job-summary info on System.err. Default value: false. This option can be set to "null" to clear the default value.

  - **REFERENCE_SEQUENCE** (str, optional) -- Reference sequence file. Default value: null.

  - **REF_FLAT** (str, required)
  - **RIBOSOMAL_INTERVALS** (str, optional)
  - **RRNA_FRAGMENT_PERCENTAGE** (int, optional)
  - **STOP_AFTER** (str, optional)
  - **STRAND_SPECIFICITY** (str, required)    - possible values: 'NONE', 'FIRST_READ_TRANSCRIPTION_STRAND', 'SECOND_READ_TRANSCRIPTION_STRAND'


  - **TMP_DIR** (str, optional) -- A file. Default value: null. This option may be specified 0 or more times.

  - **VALIDATION_STRINGENCY** (str, optional) -- Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded. Default value: STRICT. This option can be set to "null" to clear the default value.
    - possible values: 'STRICT', 'LENIENT', 'SILENT'


  - **VERBOSITY** (str, optional) -- Control verbosity of logging. Default value: INFO. This option can be set to "null" to clear the default value.
    - possible values: 'ERROR', 'WARNING', 'INFO', 'DEBUG'



**Required tools:** picard-tools

**CPU Cores:** 1

.. index:: picard_markduplicates

picard_markduplicates
=====================



    Identifies duplicate reads.

    This tool locates and tags duplicate reads in a BAM or SAM file, where
    duplicate reads are defined as originating from a single fragment of DNA.
    Duplicates can arise during sample preparation e.g. library construction
    using PCR. See also EstimateLibraryComplexity for additional notes on PCR
    duplication artifacts. Duplicate reads can also result from a single
    amplification cluster, incorrectly detected as multiple clusters by the
    optical sensor of the sequencing instrument. These duplication artifacts
    are referred to as optical duplicates.

    The MarkDuplicates tool works by comparing sequences in the 5 prime
    positions of both reads and read-pairs in a SAM/BAM file. An BARCODE_TAG
    option is available to facilitate duplicate marking using molecular
    barcodes. After duplicate reads are collected, the tool differentiates the
    primary and duplicate reads using an algorithm that ranks reads by the sums
    of their base-quality scores (default method).

    The tool's main output is a new SAM or BAM file, in which duplicates have
    been identified in the SAM flags field for each read. Duplicates are marked
    with the hexadecimal value of 0x0400, which corresponds to a decimal value
    of 1024. If you are not familiar with this type of annotation, please see
    the following blog post for additional information.

    Although the bitwise flag annotation indicates whether a read was marked as
    a duplicate, it does not identify the type of duplicate. To do this, a new
    tag called the duplicate type (DT) tag was recently added as an optional
    output in the 'optional field' section of a SAM/BAM file. Invoking the
    TAGGING_POLICY option, you can instruct the program to mark all the
    duplicates (All), only the optical duplicates (OpticalOnly), or no
    duplicates (DontTag). The records within the output of a SAM/BAM file will
    have values for the 'DT' tag (depending on the invoked TAGGING_POLICY), as
    either library/PCR-generated duplicates (LB), or sequencing-platform
    artifact duplicates (SQ). This tool uses the READ_NAME_REGEX and the
    OPTICAL_DUPLICATE_PIXEL_DISTANCE options as the primary methods to identify
    and differentiate duplicate types. Set READ_NAME_REGEX to null to skip
    optical duplicate detection, e.g. for RNA-seq or other data where duplicate
    sets are extremely large and estimating library complexity is not an aim.
    Note that without optical duplicate counts, library size estimation will be
    inaccurate.

    MarkDuplicates also produces a metrics file indicating the numbers of
    duplicates for both single- and paired-end reads.

    The program can take either coordinate-sorted or query-sorted inputs,
    however the behavior is slightly different. When the input is
    coordinate-sorted, unmapped mates of mapped records and
    supplementary/secondary alignments are not marked as duplicates. However,
    when the input is query-sorted (actually query-grouped), then unmapped
    mates and secondary/supplementary reads are not excluded from the
    duplication test and can be marked as duplicate reads.

    If desired, duplicates can be removed using the REMOVE_DUPLICATE and
    REMOVE_SEQUENCING_DUPLICATES options.

    Usage example::

        java -jar picard.jar MarkDuplicates              I=input.bam              O=marked_duplicates.bam              M=marked_dup_metrics.txt

    Documentation::

        https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates


**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/metrics**
  - **out/alignments**


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
  - **ASSUME_SORTED** (str, optional) -- If true, assume that the input file is coordinate sorted even if the header says otherwise. Default value: false. This option can be set to 'null' to clear the default value. Possible values: {true, false}
    - possible values: 'true', 'null', 'false'


  - **ASSUME_SORT_ORDER** (str, optional) -- If not null, assume that the input file has this order even if the header says otherwise. Default value: null. Possible values: {unsorted, queryname, coordinate, duplicate, unknown}
    - possible values: 'unsorted', 'queryname', 'coordinate', 'duplicate', 'unknown'


  - **BARCODE_TAG** (str, optional) -- Barcode SAM tag (ex. BC for 10X Genomics) Default value: null.

  - **COMMENT** (str, optional) -- Comment(s) to include in the output file's header. Default value: null. This option may be specified 0 or more times.

  - **COMPRESSION_LEVEL** (int, optional) -- Compression level for all compressed files created (e.g. BAM and GELI). Default value: 5.

  - **CREATE_INDEX** (bool, optional) -- Whether to create a BAM index when writing a coordinate-sorted BAM file. Default value: false. This option can be set to "null" to clear the default value. 

  - **CREATE_MD5_FILE** (bool, optional) -- Whether to create an MD5 digest for any BAM or FASTQ files created. Default value: false. This option can be set to "null" to clear the default value.

  - **DUPLICATE_SCORING_STRATEGY** (str, optional) -- The scoring strategy for choosing the non-duplicate among candidates. This option can be set to 'null' to clear the default value. Default value: SUM_OF_BASE_QUALITIES. Possible values: {SUM_OF_BASE_QUALITIES, TOTAL_MAPPED_REFERENCE_LENGTH, RANDOM}
    - possible values: 'SUM_OF_BASE_QUALITIES', 'TOTAL_MAPPED_REFERENCE_LENGTH', 'RANDOM'


  - **GA4GH_CLIENT_SECRETS** (str, optional) -- Google Genomics API client_secrets.json file path. Default value: client_secrets.json. This option can be set to "null" to clear the default value.

  - **MAX_FILE_HANDLES** (int, optional)
  - **MAX_FILE_HANDLES_FOR_READ_ENDS_MAP** (int, optional) -- Maximum number of file handles to keep open when spilling read ends to disk. Set this number a little lower than the per-process maximum number of file that may be open. This number can be found by executing the 'ulimit -n' command on a Unix system. Default value: 8000.

  - **MAX_RECORDS_IN_RAM** (int, optional) -- When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed. Default value: 500000.

  - **OPTICAL_DUPLICATE_PIXEL_DISTANCE** (int, optional) -- The maximum offset between two duplicate clusters in order to consider them optical duplicates. The default is appropriate for unpatterned versions of the Illumina platform. For the patterned flowcell models, 2500 is more appropriate. For other platforms and models, users should experiment to find what works best. Default value: 100. This option can be set to 'null' to clear the default value.

  - **PROGRAM_GROUP_COMMAND_LINE** (str, optional) -- Value of CL tag of PG record to be created. If not supplied the command line will be detected automatically. Default value: null.

  - **PROGRAM_GROUP_NAME** (str, optional) -- Value of PN tag of PG record to be created. Default value: MarkDuplicates. This option can be set to 'null' to clear the default value.

  - **PROGRAM_GROUP_VERSION** (str, optional) -- Value of VN tag of PG record to be created. If not specified, the version will be detected automatically. Default value: null.

  - **PROGRAM_RECORD_ID** (str, optional) -- The program record ID for the @PG record(s) created by this program. Set to null to disable PG record creation. This string may have a suffix appended to avoid collision with other program record IDs. Default value: MarkDuplicates. This option can be set to 'null' to clear the default value.

  - **QUIET** (bool, optional) -- Whether to suppress job-summary info on System.err. Default value: false. This option can be set to "null" to clear the default value.

  - **READ_NAME_REGEX** (str, optional) -- Regular expression that can be used to parse read names in the incoming SAM file. Read names are parsed to extract three variables: tile/region, x coordinate and y coordinate. These values are used to estimate the rate of optical duplication in order to give a more accurate estimated library size. Set this option to null to disable optical duplicate detection, e.g. for RNA-seq or other data where duplicate sets are extremely large and estimating library complexity is not an aim. Note that without optical duplicate counts, library size estimation will be inaccurate. The regular expression should contain three capture groups for the three variables, in order. It must match the entire read name. Note that if the default regex is specified, a regex match is not actually done, but instead the read name is split on colon character. For 5 element names, the 3rd, 4th and 5th elements are assumed to be tile, x and y values. For 7 element names (CASAVA 1.8), the 5th, 6th, and 7th elements are assumed to be tile, x and y values. Default value: . This option can be set to 'null' to clear the default value.

  - **READ_ONE_BARCODE_TAG** (str, optional) -- Read one barcode SAM tag (ex. BX for 10X Genomics) Default value: null.

  - **READ_TWO_BARCODE_TAG** (str, optional) -- Read two barcode SAM tag (ex. BX for 10X Genomics) Default value: null.

  - **REFERENCE_SEQUENCE** (str, optional) -- Reference sequence file. Default value: null.

  - **REMOVE_DUPLICATES** (str, optional) -- If true do not write duplicates to the output file instead of writing them with appropriate flags set. Default value: false. This option can be set to 'null' to clear the default value. Possible values: {true, false}
    - default value: true
    - possible values: 'true', 'null', 'false'


  - **REMOVE_SEQUENCING_DUPLICATES** (str, optional) -- If true remove 'optical' duplicates and other duplicates that appear to have arisen from the sequencing process instead of the library preparation process, even if REMOVE_DUPLICATES is false. If REMOVE_DUPLICATES is true, all duplicates are removed and this option is ignored. Default value: false. This option can be set to 'null' to clear the default value. Possible values: {true, false}
    - possible values: 'true', 'false', 'null'


  - **SORTING_COLLECTION_SIZE_RATIO** (float, optional) -- This number, plus the maximum RAM available to the JVM, determine the memory footprint used by some of the sorting collections. If you are running out of memory, try reducing this number. Default value: 0.25.

  - **TAGGING_POLICY** (str, optional) -- Determines how duplicate types are recorded in the DT optional attribute. Default value: DontTag. This option can be set to 'null' to clear the default value. Possible values: {DontTag, OpticalOnly, All}
    - possible values: 'DontTag', 'OpticalOnly', 'All'


  - **TAG_DUPLICATE_SET_MEMBERS** (str, optional) -- If a read appears in a duplicate set, add two tags. The first tag, DUPLICATE_SET_SIZE_TAG (DS), indicates the size of the duplicate set. The smallest possible DS value is 2 which occurs when two reads map to the same portion of the reference only one of which is marked as duplicate. The second tag, DUPLICATE_SET_INDEX_TAG (DI), represents a unique identifier for the duplicate set to which the record belongs. This identifier is the index-in-file of the representative read that was selected out of the duplicate set. Default value: false. This option can be set to 'null' to clear the default value. Possible values: {true, false}
    - possible values: 'true', 'false', 'null'


  - **TMP_DIR** (str, optional) -- A file. Default value: null. This option may be specified 0 or more times.

  - **VALIDATION_STRINGENCY** (str, optional) -- Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded. Default value: STRICT. This option can be set to "null" to clear the default value.
    - possible values: 'STRICT', 'LENIENT', 'SILENT', 'null'


  - **VERBOSITY** (str, optional) -- Control verbosity of logging. Default value: INFO. This option can be set to "null" to clear the default value.
    - possible values: 'ERROR', 'WARNING', 'INFO', 'DEBUG'



**Required tools:** picard-tools

**CPU Cores:** 12

.. index:: picard_merge_sam_bam_files

picard_merge_sam_bam_files
==========================




    Documentation::

        https://broadinstitute.github.io/picard/command-line-overview.html#MergeSamFiles


**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/alignments**


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



**Required tools:** ln (coreutils), picard-tools

**CPU Cores:** 12

.. index:: piranha

piranha
=======



    Piranha is a peak-caller for CLIP- and RIP-Seq data.
    It takes input in BED or BAM format and identifies regions of statistically significant read enrichment.
    Additional covariates may optionally be provided to further inform the peak-calling process.

    This is for the release: Piranha 1.2.0

    http://smithlabresearch.org/software/piranha/

**Input Connection**
  - **in/features**

**Output Connection**
  - **out/features**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      piranha [style=filled, fillcolor="#fce94f"];
      in_0 [label="features"];
      in_0 -> piranha;
      out_1 [label="features"];
      piranha -> out_1;
   }

**Options:**
  - **bin_size_reponse** (int, optional) -- indicates that the response (first input file) is raw reads and should be binned into bins of this size

  - **dist** (str, optional) -- Distribution type. Currently supports Poisson,NegativeBinomial, ZeroTruncatedPoisson, ZeroTruncatedNegativeBinomial (default with no covariates), PoissonRegression, NegativeBinomialRegression, ZeroTruncatedPoissonRegression, ZeroTruncatedNegativeBinomialRegression (default with covariates)
    - possible values: 'Poisson', 'NegativeBinomial', 'ZeroTruncatedPoisson', 'ZeroTruncatedNegativeBinomial', 'PoissonRegression', 'NegativeBinomialRegression', 'ZeroTruncatedPoissonRegression', 'ZeroTruncatedNegativeBinomialRegression'


  - **p_threshold** (float, optional) -- significance threshold for sites

  - **sort** (bool, optional) -- indicates that input is unsorted and Piranha should sort it for you


**Required tools:** bamtools, piranha

**CPU Cores:** 1

.. index:: pizzly

pizzly
======



    Pizzly is a tool to discover gene fusions
    in human paired-end RNA-Seq data.

    Paper:
    https://www.biorxiv.org/content/early/2017/07/20/166322

    Manual including typical usage:
    https://github.com/pmelsted/pizzly

**Input Connection**
  - **in/fusion.txt**
  - **in/second_read**
  - **in/first_read**

**Output Connection**
  - **out/unfiltered_fasta**
  - **out/fusion_fasta**
  - **out/unfiltered_json**
  - **out/log_stdout**
  - **out/log_stderr**
  - **out/fusion_json**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      pizzly [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> pizzly;
      in_1 [label="fusion.txt"];
      in_1 -> pizzly;
      in_2 [label="second_read"];
      in_2 -> pizzly;
      out_3 [label="fusion_fasta"];
      pizzly -> out_3;
      out_4 [label="fusion_json"];
      pizzly -> out_4;
      out_5 [label="log_stderr"];
      pizzly -> out_5;
      out_6 [label="log_stdout"];
      pizzly -> out_6;
      out_7 [label="unfiltered_fasta"];
      pizzly -> out_7;
      out_8 [label="unfiltered_json"];
      pizzly -> out_8;
   }

**Options:**
  - **align-score** (str, optional) -- Number of allowed mismatches

  - **cache** (str, optional) -- Cached gtf file created in a pizzly run

  - **cores** (str, required)    - default value: 6

  - **fasta-file** (str, required) -- Transcripts in .fa 

  - **gtf-file** (str, required) -- Reference transcript annotation in .gtf

  - **insert-size** (str, optional) -- Maximum insert size

  - **k-mer_size** (str, required) -- Size of k-mer


**Required tools:** pizzly

**CPU Cores:** 6

.. index:: post_cufflinksSuite

post_cufflinksSuite
===================


The cufflinks suite can be used to assembly new transcripts and
    merge those with known annotations. However, the output .gtf files
    need to be reformatted in several aspects afterwards. This step
    can be used to reformat and filter the cufflinksSuite .gtf file.

**Input Connection**
  - **in/features**

**Output Connection**
  - **out/features**
  - **out/log_stderr**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      post_cufflinksSuite [style=filled, fillcolor="#fce94f"];
      in_0 [label="features"];
      in_0 -> post_cufflinksSuite;
      out_1 [label="features"];
      post_cufflinksSuite -> out_1;
      out_2 [label="log_stderr"];
      post_cufflinksSuite -> out_2;
   }

**Options:**
  - **class-list** (str, optional) -- Class codes to be removed; possible '=,c,j,e,i,o,p,r,u,x,s,.'

  - **filter-by-class** (bool, required) -- Remove gtf if any class is found in class_code field, requieres class_list

  - **filter-by-class-and-gene-name** (bool, required) -- Combines remove-by-class and remove-by-gene-name

  - **remove-by-gene-name** (bool, required) -- Remove gtf if matches 'string' in gene_name field

  - **remove-gencode** (bool, required) -- Hard removal of gtf line which match 'ENS' in gene_name field

  - **remove-unstranded** (bool, required) -- Removes transcripts without strand specifity

  - **run_id** (str, optional) -- An arbitrary name of the new run (which is a merge of all samples).
    - default value: magic

  - **string** (str, optional) -- String to match in gtf field gene_name for discarding


**Required tools:** cat (coreutils)

**CPU Cores:** 6

.. index:: post_sawdust

post_sawdust
============



    bla bla

**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/alignments**
  - **out/log_stderr**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      post_sawdust [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> post_sawdust;
      out_1 [label="alignments"];
      post_sawdust -> out_1;
      out_2 [label="log_stderr"];
      post_sawdust -> out_2;
   }

**Options:**
  - **library_type** (str, required)    - possible values: 'fr-unstranded', 'fr-firststrand', 'fr-secondstrand'


  - **prf** (bool, optional)
  - **read_type** (str, required)    - possible values: 'single', 'paired'


  - **seq_type** (str, required)    - possible values: 'RNA', 'DNA'


  - **split_ident** (str, required)    - default value:  


**Required tools:** cat (coreutils), pigz, samtools

**CPU Cores:** 2

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
    file, then simply type::

        preseq c_curve -o output.txt input.sort.bed

    Documentation::

        http://smithlabresearch.org/software/preseq/


**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/complexity_curve**


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

  - **step** (int, optional) -- step size in extrapolations (default: 1e+06)

  - **vals** (bool, optional) -- input is a text file containing only the observed counts

  - **verbose** (bool, optional) -- print more information


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

**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/future_genome_coverage**


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

**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/future_yield**


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

.. index:: reformatCigar

reformatCigar
=============



    The segemehl release from April 2017 (segemehl/export.13) uses the
    expanded version of Cigar strings: It clearly separates between match (=)
    and mismatch (X) while in older Cigar string definitions both are
    encoded with an M character.
    However, subsequent tools like htseq-count are not able to handle these
    newer versions of Cigar strings.

    This step transforms the new Cigar string encoding back to the old one.

    This step depends on the segemehl_2017 step. You still may need to
    add the subsequenct step s2c before calling cufflinks and/or htseq-count.

**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/alignments**
  - **out/log**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      reformatCigar [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> reformatCigar;
      out_1 [label="alignments"];
      reformatCigar -> out_1;
      out_2 [label="log"];
      reformatCigar -> out_2;
   }

**Options:**
  - **blocksize** (int, optional) -- Blocksize to read the input file, in Megabytes.Default: 2 (2,000,000 bytes)

  - **threads** (int, optional) -- Number of threads 2B started. (Default: 1). Beware that this is only for (un-)compressing, the reformating is using a single CPU only.


**Required tools:** cat (coreutils), pigz

**CPU Cores:** 1

.. index:: remove_duplicate_reads_runs

remove_duplicate_reads_runs
===========================



    Duplicates are removed by Picard tools 'MarkDuplicates'.

    typical command line::

        MarkDuplicates INPUT=<SAM/BAM> OUTPUT=<SAM/BAM>
                       METRICS_FILE=<metrics-out> REMOVE_DUPLICATES=true

**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/metrics**
  - **out/alignments**


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

**Required tools:** MarkDuplicates

**CPU Cores:** 12

.. index:: rgt_thor

rgt_thor
========



    THOR is an HMM-based approach to detect and analyze differential peaks in
    two sets of ChIP-seq data from distinct biological conditions with
    replicates. THOR performs genomic signal processing, peak calling and
    p-value calculation in an integrated framework. For differential peak
    calling without replicates use ODIN.

    For more information please refer to:

    Allhoff, M., Sere K., Freitas, J., Zenke, M.,  Costa, I.G. (2016),
    Differential Peak Calling of ChIP-seq Signals with Replicates with THOR,
    Nucleic Acids Research, epub gkw680 [paper][supp].

    Feel free to post your question in our googleGroup or write an e-mail:
    rgtusers@googlegroups.com

**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/diff_peaks_bed**
  - **out/thor_config**
  - **out/thor_setup_info**
  - **out/diff_narrow_peaks**
  - **out/chip_seq_bigwig**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      rgt_thor [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> rgt_thor;
      out_1 [label="chip_seq_bigwig"];
      rgt_thor -> out_1;
      out_2 [label="diff_narrow_peaks"];
      rgt_thor -> out_2;
      out_3 [label="diff_peaks_bed"];
      rgt_thor -> out_3;
      out_4 [label="thor_config"];
      rgt_thor -> out_4;
      out_5 [label="thor_setup_info"];
      rgt_thor -> out_5;
   }

**Options:**
  - **binsize** (int, optional) -- Size of bins for creating the signal.

  - **chrom_sizes_file** (str, required)
  - **config_file** (dict, required) -- A dictionary with 

  - **deadzones** (str, optional) -- Define blacklisted genomic regions to be ignored by the peak caller. 

  - **exts** (str, optional) -- Read's extension size for BAM files (comma separated list for each BAM file in config file). If option is not chosen, estimate extension sizes from reads.

  - **factors-inputs** (str, optional) -- Normalization factors for input-DNA (comma separated list for each BAM file in config file). If option is not chosen, estimate factors.

  - **genome** (str, required) -- FASTA file containing the complete genome sequence

  - **housekeeping-genes** (str, optional) -- Define housekeeping genes (BED format) used for normalizing.

  - **merge** (bool, optional) -- Merge peaks which have a distance less than the estimated mean fragment size (recommended for histone data).

  - **no-correction** (bool, optional) -- Do not use multiple test correction for p-values (Benjamini/Hochberg).

  - **no-gc-content** (bool, optional) -- Do not normalize towards GC content.

  - **pvalue** (float, optional) -- P-value cutoff for peak detection. Call only peaks with p-value lower than cutoff.

  - **report** (bool, optional) -- Generate HTML report about experiment.

  - **save-input** (bool, optional) -- Save input DNA bigwig (if input was provided).

  - **scaling-factors** (str, optional) -- Scaling factor for each BAM file (not control input-DNA) as comma separated list for each BAM file in config file. If option is not chosen, follow normalization strategy (TMM or HK approach)

  - **step** (int, optional) -- Stepsize with which the window consecutively slides across the genome to create the signal.


**Required tools:** printf (coreutils), rgt-THOR

**CPU Cores:** 4

.. index:: rseqc

rseqc
=====



    The RSeQC step can be used to evaluate aligned reads in a BAM file. RSeQC
    does not only report raw sequence-based metrics, but also quality control
    metrics like read distribution, gene coverage, and sequencing depth.

**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/DupRate_stdout**
  - **out/geneBody_coverage_stdout**
  - **out/bam_stat**
  - **out/junctionSaturation_r**
  - **out/DupRate_seq**
  - **out/inner_distance** (optional)
  - **out/junction_saturation_stderr**
  - **out/gc_r**
  - **out/gc_stderr**
  - **out/geneBody_coverage.r**
  - **out/read_distribution**
  - **out/junction_plot**
  - **out/inner_distance_stderr** (optional)
  - **out/gc_xls**
  - **out/inner_distance_stdout** (optional)
  - **out/inner_distance_freq** (optional)
  - **out/geneBody_coverage.txt**
  - **out/geneBody_coverage_stderr**
  - **out/junction_bed**
  - **out/junction_annotation_stdout**
  - **out/gc_stdout**
  - **out/infer_experiment**
  - **out/inner_distance_plot** (optional)
  - **out/DupRate_plot_r**
  - **out/DupRate_pos**
  - **out/DupRate_stderr**
  - **out/junction_xls**
  - **out/junction_saturation_stdout**
  - **out/junction_annotation_stderr**


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
      out_1 [label="DupRate_plot_r"];
      rseqc -> out_1;
      out_2 [label="DupRate_pos"];
      rseqc -> out_2;
      out_3 [label="DupRate_seq"];
      rseqc -> out_3;
      out_4 [label="DupRate_stderr"];
      rseqc -> out_4;
      out_5 [label="DupRate_stdout"];
      rseqc -> out_5;
      out_6 [label="bam_stat"];
      rseqc -> out_6;
      out_7 [label="gc_r"];
      rseqc -> out_7;
      out_8 [label="gc_stderr"];
      rseqc -> out_8;
      out_9 [label="gc_stdout"];
      rseqc -> out_9;
      out_10 [label="gc_xls"];
      rseqc -> out_10;
      out_11 [label="geneBody_coverage.r"];
      rseqc -> out_11;
      out_12 [label="geneBody_coverage.txt"];
      rseqc -> out_12;
      out_13 [label="geneBody_coverage_stderr"];
      rseqc -> out_13;
      out_14 [label="geneBody_coverage_stdout"];
      rseqc -> out_14;
      out_15 [label="infer_experiment"];
      rseqc -> out_15;
      out_16 [label="inner_distance", style=filled, fillcolor="#a7a7a7"];
      rseqc -> out_16;
      out_17 [label="inner_distance_freq", style=filled, fillcolor="#a7a7a7"];
      rseqc -> out_17;
      out_18 [label="inner_distance_plot", style=filled, fillcolor="#a7a7a7"];
      rseqc -> out_18;
      out_19 [label="inner_distance_stderr", style=filled, fillcolor="#a7a7a7"];
      rseqc -> out_19;
      out_20 [label="inner_distance_stdout", style=filled, fillcolor="#a7a7a7"];
      rseqc -> out_20;
      out_21 [label="junctionSaturation_r"];
      rseqc -> out_21;
      out_22 [label="junction_annotation_stderr"];
      rseqc -> out_22;
      out_23 [label="junction_annotation_stdout"];
      rseqc -> out_23;
      out_24 [label="junction_bed"];
      rseqc -> out_24;
      out_25 [label="junction_plot"];
      rseqc -> out_25;
      out_26 [label="junction_saturation_stderr"];
      rseqc -> out_26;
      out_27 [label="junction_saturation_stdout"];
      rseqc -> out_27;
      out_28 [label="junction_xls"];
      rseqc -> out_28;
      out_29 [label="read_distribution"];
      rseqc -> out_29;
   }

**Options:**
  - **reference** (str, required) -- Reference gene model in bed fomat. [required]

  - **treatAs** (str, required) -- Some modules in rseqc need paired end dataan fail otherwise on single end [required]
    - possible values: 'single', 'paired'



**Required tools:** bam_stat.py, cat (coreutils), geneBody_coverage.py, infer_experiment.py, inner_distance.py, junction_annotation.py, junction_saturation.py, read_GC.py, read_distribution.py, read_duplication.py, rm (coreutils)

**CPU Cores:** 4

.. index:: s2c

s2c
===



    s2c formats the output of segemehl mapping to be compatible with
    the cufflinks suite of tools for differential expr. analysis of
    RNA-Seq data and their visualisation.
    For details on cufflinks we refer to the author's webpage:

    http://cole-trapnell-lab.github.io/cufflinks/


**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/alignments**
  - **out/log**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      s2c [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> s2c;
      out_1 [label="alignments"];
      s2c -> out_1;
      out_2 [label="log"];
      s2c -> out_2;
   }

**Options:**
  - **maxDist** (int, optional) -- specifies the maximal distance of a splice junction. junctions with disctance higher than this value are classified as fusions (default is 200.000nt)

  - **tmp_dir** (str, required) -- Temp directory for 's2c.py'. This can be in the /work/username/ path, since it is only temporary.


**Required tools:** cat (coreutils), dd (coreutils), pigz, samtools

**CPU Cores:** 6

.. index:: salmon

salmon
======



    Salmon



**Input Connection**
  - **in/second_read** (optional)
  - **in/first_read**

**Output Connection**
  - **out/lib_format_counts.json**
  - **out/quant.genes.sf**
  - **out/flenDist.txt**
  - **out/ambig_info.tsv**
  - **out/observed_bias_3p.gz**
  - **out/meta_info.json**
  - **out/observed_bias.gz**
  - **out/fld.gz**
  - **out/log_stdout**
  - **out/log_stderr**
  - **out/quant.sf**
  - **out/expected_bias.gz**
  - **out/cmd_info.json**
  - **out/salmon_quant.log**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      salmon [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> salmon;
      in_1 [label="second_read", style=filled, fillcolor="#a7a7a7"];
      in_1 -> salmon;
      out_2 [label="ambig_info.tsv"];
      salmon -> out_2;
      out_3 [label="cmd_info.json"];
      salmon -> out_3;
      out_4 [label="expected_bias.gz"];
      salmon -> out_4;
      out_5 [label="fld.gz"];
      salmon -> out_5;
      out_6 [label="flenDist.txt"];
      salmon -> out_6;
      out_7 [label="lib_format_counts.json"];
      salmon -> out_7;
      out_8 [label="log_stderr"];
      salmon -> out_8;
      out_9 [label="log_stdout"];
      salmon -> out_9;
      out_10 [label="meta_info.json"];
      salmon -> out_10;
      out_11 [label="observed_bias.gz"];
      salmon -> out_11;
      out_12 [label="observed_bias_3p.gz"];
      salmon -> out_12;
      out_13 [label="quant.genes.sf"];
      salmon -> out_13;
      out_14 [label="quant.sf"];
      salmon -> out_14;
      out_15 [label="salmon_quant.log"];
      salmon -> out_15;
   }

**Options:**
  - **cores** (int, optional) -- workaround to specify cores for grid engine and threads ie
    - default value: 1

  - **g** (str, optional) -- File containing a mapping of transcripts to genes. If this file is provided Salmon will output both quant.sf and quant.genes.sf files, where the latter contains aggregated gene-level abundance estimates. The transcript to gene mapping should be provided as either a GTF file, or a in a simple tab-delimited format where each line contains the name of a transcript and the gene to which it belongs separated by a tab. The extension of the file is used to determine how the file should be parsed. Files ending in '.gtf', '.gff' or '.gff3'are assumed to be in GTF format; files with any other extension are assumed to be in the simple format. In GTF / GFF format, the 'transcript_id' is assumed to contain the transcript identifier and the 'gene_id' is assumed to contain the corresponding gene identifier.

  - **i** (str, required) -- Salmon index


**Required tools:** mkdir (coreutils), mv (coreutils), rm (coreutils), salmon

**CPU Cores:** 1

.. index:: sam_to_fastq

sam_to_fastq
============



    Tbla bla bla    http://www.htslib.org/doc/samtools.html

**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/first_read**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      sam_to_fastq [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> sam_to_fastq;
      out_1 [label="first_read"];
      sam_to_fastq -> out_1;
   }

**Options:**
  - **F** (int, optional)
  - **addF** (int, optional)
  - **f** (int, optional)

**Required tools:** pigz, samtools

**CPU Cores:** 8

.. index:: sam_to_sorted_bam

sam_to_sorted_bam
=================



    The step sam_to_sorted_bam builds on 'samtools sort' to sort SAM files and
    output BAM files.

    Sort alignments by leftmost coordinates, or by read name when -n is used.
    An appropriate @HD-SO sort order header tag will be added or an existing
    one updated if necessary.

    Documentation::

        http://www.htslib.org/doc/samtools.html

**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/alignments**


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
  - **dd-blocksize** (str, optional)    - default value: 256k

  - **genome-faidx** (str, required)
  - **sort-by-name** (bool, required)
  - **temp-sort-dir** (str, optional) -- Intermediate sort files are stored intothis directory.


**Required tools:** dd (coreutils), pigz, samtools

**CPU Cores:** 8

.. index:: samtools

samtools
========



    The step samtools wraps parts of the 'samtools' packages. It is intended for
    reformatting SAM/BAM files and not completely implemented.

    Feel free to add the samtools options you need!

    The options listed below are implemented.

    For a description/explanation about SAM flags, we refer to
    - the samtools manual page http://www.htslib.org/doc/samtools.html
    - the Picard page to explain SAM flags https://broadinstitute.github.io/picard/explain-flags.html


**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/alignments**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      samtools [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> samtools;
      out_1 [label="alignments"];
      samtools -> out_1;
   }

**Options:**
  - **F_skip** (int, optional) -- Do not output alignments with any bits set in INT present in the FLAG field.

  - **dd-blocksize** (str, optional) -- Blocksize for dd tool.
    - default value: 2048

  - **f_keep** (int, optional) -- Only output alignments with all bits set in INT present in the FLAG field.

  - **genome-faidx** (str, required)
  - **keep_header** (bool, optional) -- Include the header in the output.
    - default value: True

  - **output_bam** (bool, optional) -- Output in the BAM format.

  - **q_mapq** (int, optional) -- Skip alignments with MAPQ smaller than this value.

  - **sort-by-name** (bool, required) -- Sort by read names (i.e., the QNAME field) rather than by chromosomal coordinates.

  - **temp-sort-dir** (str, required) -- Intermediate sort files are stored intothis directory.


**Required tools:** dd (coreutils), pigz, samtools

**CPU Cores:** 8

.. index:: samtools_faidx

samtools_faidx
==============



    Index reference sequence in the FASTA format or extract subsequence from
    indexed reference sequence. If no region is specified, faidx will index the
    file and create <ref.fasta>.fai on the disk. If regions are specified, the
    subsequences will be retrieved and printed to stdout in the FASTA format.

    The sequences in the input file should all have different names. If they do
    not, indexing will emit a warning about duplicate sequences and retrieval
    will only produce subsequences from the first sequence with the duplicated
    name.

**Input Connection**
  - **in/sequence**

**Output Connection**
  - **out/indices**


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

**Required tools:** mv (coreutils), samtools

**CPU Cores:** 4

.. index:: samtools_index

samtools_index
==============



    Index a coordinate-sorted BAM or CRAM file for fast random access.
    (Note that this does not work with SAM files even if they are bgzip
    compressed to index such files, use tabix(1) instead.)

    Documentation::

        http://www.htslib.org/doc/samtools.html

**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/alignments**
  - **out/indices**
  - **out/index_stats**


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
  - **index_type** (str, required)    - possible values: 'bai', 'csi'



**Required tools:** ln (coreutils), samtools

**CPU Cores:** 4

.. index:: samtools_merge

samtools_merge
==============


The step samtools_merge builds on 'samtools merge' to merge sorted SAM/BAM files and output
    SAM/BAM/CRAM files.

    Merge adds readgroup tag to preserve the information on the orginial sample.

    This step wraps samtools merge from release 1.7.

    Documentation:

        http://www.htslib.org/doc/samtools.html


**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/alignments**
  - **out/log**
  - **out/err**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      samtools_merge [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> samtools_merge;
      out_1 [label="alignments"];
      samtools_merge -> out_1;
      out_2 [label="err"];
      samtools_merge -> out_2;
      out_3 [label="log"];
      samtools_merge -> out_3;
   }

**Options:**
  - **1** (bool, optional) -- compress level 1

  - **R** (str, optional) -- merge file in the specified region STR [all]

  - **c** (bool, optional) -- Combine @RG headers with colliding IDs [alter IDs to be distinct]

  - **dd-blocksize** (str, optional)    - default value: 4M

  - **f** (bool, optional) -- overwrite the output BAM if exist

  - **l** (int, optional) -- compression level, from 0 to 9 [-1]

  - **n** (bool, optional) -- sort by read names

  - **p** (bool, optional) -- Combine @PG headers with colliding IDs [alter IDs to be distinct]

  - **pigz-blocksize** (str, optional)    - default value: 4096

  - **r** (bool, optional) -- attach RG tag (inferred from file names)

  - **run_id** (str, optional) -- A name for the run. Since this step merges multiple samples into a single one, the run_id cannot be the sample name anymore.
    - default value: mergeSBam

  - **s** (str, optional) -- override random seed

  - **t** (str, optional) -- Input files are sorted by TAG value

  - **threads** (int, optional) -- Number of additional threads to use [0]

  - **u** (bool, optional) -- uncompressed BAM output


**Required tools:** dd (coreutils), pigz, samtools

**CPU Cores:** 6

.. index:: samtools_sort

samtools_sort
=============



    automatically recognizes input format
    The step implements samtools sort to sort sam, cram, bam.
    crashes on sam to cram

    Sort alignments by leftmost coordinates, or by read name when -n is used.
    An appropriate @HD-SO sort order header tag will be added or an existing
    one updated if necessary.

    Documentation::

        http://www.htslib.org/doc/samtools.html

**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/alignments**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      samtools_sort [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> samtools_sort;
      out_1 [label="alignments"];
      samtools_sort -> out_1;
   }

**Options:**
  - **O** (str, required) -- output format bam sam or cram
    - possible values: 'BAM', 'SAM', 'CRAM'


  - **cores** (int, optional) -- workaround to specify cores for grid engine and threads ie
    - default value: 1

  - **dd-blocksize** (str, optional) -- Read data with ``dd`` and set the blocksize.
    - default value: 4096k

  - **fifo** (bool, optional) -- Enable the FIFO functionality for splitting large input files.

  - **l** (str, optional) -- Set compression level, from 0 (uncompressed) to 9 (best)

  - **m** (str, optional) -- Set maximum memory per thread; suffix K/M/G recognized [768M]

  - **reference** (str, optional) -- reference fasta file need for cram output

  - **sort-by-name** (bool, required) -- sort by read name original option is -n

  - **t** (str, optional) -- Sort by value of TAG. Uses position as secondary index (or read name if -n is set)

  - **temp-sort-dir** (str, optional) -- Intermediate sort files are stored intothis directory. original option -T


**Required tools:** dd (coreutils), pigz, samtools

**CPU Cores:** 1

.. index:: samtools_stats

samtools_stats
==============



    samtools stats collects statistics from BAM files and outputs in a text
    format. The output can be visualized graphically using plot-bamstats.

    Documentation::

        http://www.htslib.org/doc/samtools.html

**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/stats**


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
  - **dd-blocksize** (str, optional)    - default value: 256k


**Required tools:** dd (coreutils), pigz, samtools

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
    unmapped reads::

       mkfifo genome_fifo unmapped_fifo
       cat <genome-fasta> -o genome_fifo

    The executed segemehl command is this::

        segemehl -d genome_fifo -i <genome-index-file> -q <read1-fastq>
                 [-p <read2-fastq>] -u unmapped_fifo -H 1 -t 11 -s -S -D 0
                 -o /dev/stdout |  pigz --blocksize 4096 --processes 2 -c

    The unmapped reads are saved via these commands::

        cat unmapped_fifo | pigz --blocksize 4096 --processes 2 -c >
        <unmapped-fastq>


**Input Connection**
  - **in/second_read**
  - **in/first_read**

**Output Connection**
  - **out/alignments**
  - **out/log**
  - **out/unmapped**


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

  - **dd-blocksize** (str, optional)    - default value: 2M

  - **differences** (int, optional) -- search seeds initially with <n> differences (default:1)
    - default value: 1

  - **dropoff** (int, optional) -- dropoff parameter for extension (default:8)

  - **evalue** (float, optional) -- max evalue (default:5.000000)

  - **extensionpenalty** (int, optional) -- penalty for a mismatch during extension (default:4)

  - **extensionscore** (int, optional) -- score of a match during extension (default:2)

  - **filebins** (str, optional) -- file bins with basename <string> for easier data handling (default:none)
    - default value: none

  - **fix-qnames** (bool, optional) -- The QNAMES field of the input will be purged from spaces and everything thereafter.

  - **genome** (str, required) -- Path to genome file

  - **hardclip** (bool, optional) -- enable hard clipping

  - **hitstrategy** (int, optional) -- report only best scoring hits (=1) or all (=0) (default:1)
    - default value: 1
    - possible values: '0', '1'


  - **index** (str, required) -- path/filename of db index (default:none)

  - **index2** (str, optional) -- path/filename of second db index (default:none)

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

  - **pigz-blocksize** (str, optional)    - default value: 2048

  - **polyA** (bool, optional) -- clip polyA tail

  - **prime3** (str, optional) -- add 3' adapter (default:none)

  - **prime5** (str, optional) -- add 5' adapter (default:none)

  - **showalign** (bool, optional) -- show alignments

  - **silent** (bool, optional) -- shut up!
    - default value: True

  - **splicescorescale** (float, optional) -- report spliced alignment with score s only if <f>*s is larger than next best spliced alignment (default:1.000000)

  - **splits** (bool, optional) -- detect split/spliced reads (default:none)

  - **threads** (int, optional) -- start <n> threads (default:10)
    - default value: 10


**Required tools:** cat (coreutils), dd (coreutils), mkfifo (coreutils), pigz, segemehl

**CPU Cores:** 10

.. index:: segemehl_2017

segemehl_2017
=============



    segemehl_2017 is a software to map short sequencer reads to reference genomes.
    Unlike other methods, segemehl is able to detect not only mismatches but
    also insertions and deletions. Furthermore, segemehl is not limited to a
    specific read length and is able to mapprimer- or polyadenylation
    contaminated reads correctly.

    This step is a wrapper for an unpublished version of segemehl that we re-
    ceived from Steve Hoffmann on April, 24th 2017.
    It automatically generates the bed file containing the realigned split
    reads (split junctions).
    No need to run testrealign.x after segemehl anymore.

    This step creates at first two FIFOs. The first is used to provide the
    genome data for segemehl and the second is used for the output of the
    unmapped reads::

       mkfifo genome_fifo unmapped_fifo
       cat <genome-fasta> -o genome_fifo

    The executed segemehl command is this::

        segemehl -d genome_fifo -i <genome-index-file> -q <read1-fastq>
                 [-p <read2-fastq>] -u unmapped_fifo -H 1 -t 11 -s -S -D 0
                 -o /dev/stdout |  pigz --blocksize 4096 --processes 2 -c

    The unmapped reads are saved via these commands::

        cat unmapped_fifo | pigz --blocksize 4096 --processes 2 -c >
        <unmapped-fastq>


**Input Connection**
  - **in/second_read**
  - **in/first_read**

**Output Connection**
  - **out/splits_sngl**
  - **out/unmapped**
  - **out/alignments**
  - **out/splits_trns**
  - **out/log**
  - **out/splits_mult**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      segemehl_2017 [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> segemehl_2017;
      in_1 [label="second_read"];
      in_1 -> segemehl_2017;
      out_2 [label="alignments"];
      segemehl_2017 -> out_2;
      out_3 [label="log"];
      segemehl_2017 -> out_3;
      out_4 [label="splits_mult"];
      segemehl_2017 -> out_4;
      out_5 [label="splits_sngl"];
      segemehl_2017 -> out_5;
      out_6 [label="splits_trns"];
      segemehl_2017 -> out_6;
      out_7 [label="unmapped"];
      segemehl_2017 -> out_7;
   }

**Options:**
  - **MEOP** (bool, optional) -- output MEOP field for easier variance calling in SAM (XE:Z:)

  - **accuracy** (int, optional) -- min percentage of matches per read in semi-global alignment (default:90)

  - **bisulfite** (int, optional) -- bisulfite mapping with methylC-seq/Lister et al. (=1) or bs-seq/Cokus et al. protocol (=2) (default:0)
    - possible values: '0', '1', '2'


  - **brief** (bool, optional) -- brief output

  - **briefcigar** (bool, optional) -- brief cigar string (M vs X and =)

  - **checkidx** (bool, optional) -- check index

  - **clipacc** (int, optional) -- clipping accuracy (default:70)

  - **database** (str, required) -- (Space separated list of ) filename(s) of database (e.g. genome) sequene(s)

  - **dd-blocksize** (str, optional)    - default value: 1M

  - **differences** (int, optional) -- search seeds initially with <n> differences (default:1)

  - **dropoff** (int, optional) -- dropoff parameter for extension (default:8)

  - **evalue** (float, optional) -- max evalue (default:5.000000)

  - **extensionpenalty** (int, optional) -- penalty for a mismatch during extension (default:4)

  - **filebins** (str, optional) -- file bins with basename <string> for easier data handling (default:none)

  - **fix-qnames** (bool, optional) -- The QNAMES field of the input will be purged from spaces and everything thereafter.

  - **hitstrategy** (int, optional) -- report only best scoring hits (=1) or all (=0) (default:1)
    - possible values: '0', '1'


  - **index** (str, optional) -- Path to database index for segemehl (default:none)

  - **index2** (str, optional) -- Path to second database index for segemehl (default:none)

  - **jump** (int, optional) -- search seeds with jump size <n> (0=automatic) (default:0)

  - **maxinsertsize** (int, optional) -- maximum size of the inserts (paired end) (default:5000)

  - **maxinterval** (int, optional) -- maximum width of a suffix array interval, i.e. a query seed will be omitted if it matches more than <n> times (default:100)

  - **maxout** (int, optional) -- maximum number of alignments that will be reported. If set to zero, all alignments will be reported (default:0)

  - **maxsplitevalue** (float, optional) -- max evalue for splits (default:50.000000)

  - **minfraglen** (int, optional) -- min length of a spliced fragment (default:20)

  - **minfragscore** (int, optional) -- min score of a spliced fragment (default:18)

  - **minsize** (int, optional) -- minimum size of queries (default:12)

  - **minsplicecover** (int, optional) -- min coverage for spliced transcripts (default:80)

  - **nohead** (bool, optional) -- do not output header

  - **nosuflinks** (bool, optional) -- dont use suflinks (does not affect index construction, for short reads only, increases runtime!)

  - **prime3** (str, optional) -- add 3' adapter (default:none)

  - **prime5** (str, optional) -- add 5' adapter (default:none)

  - **readgroupfile** (str, optional) -- filename to read @RG header (default:none)

  - **readgroupid** (str, optional) -- read group id (default:none)

  - **showalign** (bool, optional) -- show alignments

  - **splicescorescale** (float, optional) -- report spliced alignment with score s only if <f>*s is larger than next best spliced alignment (default:1.000000)

  - **splits** (bool, optional) -- detect split/spliced reads (default:none)

  - **threads** (int, optional) -- start <n> threads (default:1)
    - default value: 1


**Required tools:** cat (coreutils), dd (coreutils), mkfifo (coreutils), pigz, segemehl

**CPU Cores:** 10

.. index:: segemehl_generate_index

segemehl_generate_index
=======================



    The step segemehl_generate_index generates a index for given reference
    sequences.

    Documentation::

       http://www.bioinf.uni-leipzig.de/Software/segemehl/

**Input Connection**
  - **in/reference_sequence**

**Output Connection**
  - **out/log**
  - **out/segemehl_index**


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
  - **dd-blocksize** (str, optional)    - default value: 2M

  - **index-basename** (str, required) -- Basename for created segemehl index.

  - **pigz-blocksize** (str, optional)    - default value: 2048

  - **threads** (int, optional) -- start <n> threads (default:4)


**Required tools:** dd (coreutils), mkfifo (coreutils), pigz, segemehl

**CPU Cores:** 4

.. index:: segemehl_generate_index_bisulfite

segemehl_generate_index_bisulfite
=================================



    The step segemehl_generate_index_bisulfite generates a pair of
    indexes for given reference for segemehl bisulfite sequencing
    mapping.

    Documentation::

       http://www.bioinf.uni-leipzig.de/Software/segemehl/

**Input Connection**
  - **in/reference_sequence**

**Output Connection**
  - **out/indexCT**
  - **out/log**
  - **out/indexGA**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      segemehl_generate_index_bisulfite [style=filled, fillcolor="#fce94f"];
      in_0 [label="reference_sequence"];
      in_0 -> segemehl_generate_index_bisulfite;
      out_1 [label="indexCT"];
      segemehl_generate_index_bisulfite -> out_1;
      out_2 [label="indexGA"];
      segemehl_generate_index_bisulfite -> out_2;
      out_3 [label="log"];
      segemehl_generate_index_bisulfite -> out_3;
   }

**Options:**
  - **bisulfite** (int, optional) -- bisulfite mapping with methylC-seq/Lister et al. (=1) or bs-seq/Cokus et al. protocol (=2) (default:0)
    - default value: 1
    - possible values: '1', '2'


  - **dd-blocksize** (str, optional)    - default value: 2M

  - **generate** (str, optional) -- Filename of first (CT) db index that is generated and store to disk (efault: index-basename.ctidx).

  - **generate2** (str, optional) -- Filename of 2nd (GA) db index that is generated and store to disk (efault: index-basename.ctidx).

  - **index-basename** (str, required) -- Basename for created segemehl index.

  - **pigz-blocksize** (str, optional)    - default value: 2048

  - **threads** (int, optional) -- start <n> threads (default:4)


**Required tools:** dd (coreutils), mkfifo (coreutils), pigz, segemehl

**CPU Cores:** 4

.. index:: soapfuse

soapfuse
========



    SOAPfuse is a tool to discover gene fusions
    in human paired-end RNA-Seq data.

    Paper:
    https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-2-read12

    Manual including required folder structure and typical usage:
    https://sourceforge.net/p/soapfuse/wiki/Home/

**Input Connection**
  - **in/second_read**
  - **in/first_read**

**Output Connection**
  - **out/tar_archive**
  - **out/log_stdout**
  - **out/sf_sample_list**
  - **out/log_stderr**
  - **out/sf_config**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      soapfuse [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> soapfuse;
      in_1 [label="second_read"];
      in_1 -> soapfuse;
      out_2 [label="log_stderr"];
      soapfuse -> out_2;
      out_3 [label="log_stdout"];
      soapfuse -> out_3;
      out_4 [label="sf_config"];
      soapfuse -> out_4;
      out_5 [label="sf_sample_list"];
      soapfuse -> out_5;
      out_6 [label="tar_archive"];
      soapfuse -> out_6;
   }

**Options:**
  - **c** (str, required) -- SOAPfuse config; In the config file following variables are overwritten: path to index: DB_db_dir path to soapfuse bin: PG_pg_dir path to soapfuse source: PS_ps_dir suffix for fastq: PA_all_fq_postfix (i.e.: \*fastq.gz) cores: PA_all_process_of_align_software 

  - **cores** (int, required)    - default value: 6

  - **es** (int, optional) -- The step you want to end at 1-9
    - default value: 8

  - **path_to_index_dir** (str, required) -- Sets 'DB_db_dir' in SOAPfuse config

  - **path_to_sf_bin_dir** (str, required) -- Sets 'PG_pg_dir' in SOAPfuse config

  - **path_to_sf_source** (str, required) -- Sets 'PS_ps_dir' in SOAPfuse config

  - **read_length** (int, required) -- Sets read length for the sample list

  - **suffix_for_fq_file** (str, required) -- Sets 'PA_all_fq_postfix' in SOAPfuse config


**Required tools:** cp (coreutils), echo, ln (coreutils), mkdir (coreutils), rm (coreutils), soapfuse, tar

**CPU Cores:** 6

.. index:: source_controller

source_controller
=================



    This step combines all inputs, produces a symlink to each file
    and hashes them. It may be use to inherit from source steps
    so changes in the source files can be detected later on.

**Input Connection**
  - **in/raw** - Files to control.

**Output Connection**
  - **out/merged** - All controlled files combined in one run ``links``. The output files are named ``<previous run id>-<file name>``.


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      source_controller [style=filled, fillcolor="#fce94f"];
      in_0 [label="raw"];
      in_0 -> source_controller;
      out_1 [label="merged"];
      source_controller -> out_1;
   }

**Options:**
  - **cores** (int, optional) -- Number of threads used to calculate the hash sums.
    - default value: 4


**Required tools:** ln (coreutils)

**CPU Cores:** 4

.. index:: split_fastq

split_fastq
===========





**Input Connection**
  - **in/second_read**
  - **in/first_read**

**Output Connection**
  - **out/first_read**
  - **out/second_read**
  - **out/log_stdout**
  - **out/log_stderr**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      split_fastq [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> split_fastq;
      in_1 [label="second_read"];
      in_1 -> split_fastq;
      out_2 [label="first_read"];
      split_fastq -> out_2;
      out_3 [label="log_stderr"];
      split_fastq -> out_3;
      out_4 [label="log_stdout"];
      split_fastq -> out_4;
      out_5 [label="second_read"];
      split_fastq -> out_5;
   }

**Options:**
  - **cores** (int, optional) -- workaround to specify cores for grid engine and threads ie
    - default value: 1

  - **outfile_count** (int, required) -- Number of outfiles

  - **readcount** (int, required) -- Number of reads per targetfile


**Required tools:** split_fastqn

**CPU Cores:** 1

.. index:: sra_fastq_dump

sra_fastq_dump
==============



    sra tools is a suite from NCBI to handle sra (short read archive) files.
    fastq-dump is an sra tool that dumps the content of an sra file in fastq
    format

    The following options cannot be set, as they would interefere with the
    pipeline implemented in this step::

        -O|--outdir <path>
            Output directory, default is working directory '.'
        -Z|--stdout
            Output to stdout, all split data become joined into single stream
        --gzip
            Compress output using gzip
        --bzip2
            Compress output using bzip2

    **Multiple File Options**
    Setting these options will produce more
    than 1 file, each of which will be suffixed
    according to splitting criteria::

        --split-files
            Dump each read into separate file.Files
            will receive suffix corresponding to read number
        --split-3
            Legacy 3-file splitting for mate-pairs:
            First biological reads satisfying dumping
            conditions are placed in files \*_1.fastq and
            \*_2.fastq If only one biological read is
            present it is placed in \*.fastq Biological
            reads and above are ignored.
        -G|--spot-group
            Split into files by SPOT_GROUP (member name)
        -R|--read-filter <[filter]>
            Split into files by READ_FILTER value
            optionally filter by value:
            pass|reject|criteria|redacted
        -T|--group-in-dirs
            Split into subdirectories instead of files
        -K|--keep-empty-files
            Do not delete empty files


    Details on fastq-dump can be found at
    https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump

    To make IO cluster friendly, fastq-dump is not reading th sra file directly.
    Rather, dd with configurable blocksize is used to provide the sra file via a
    fifo to fastq-dump.

    The executed calls lools like this

    mkfifo sra_fifo
    dd bs=4M if=<sra-file> of=sra_fifo
    fastq-dump -Z sra_fifo | pigz --blocksize 4096 --processes 2 > file.fastq

**Input Connection**
  - **in/sequence**

**Output Connection**
  - **out/second_read**
  - **out/log**
  - **out/first_read**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      sra_fastq_dump [style=filled, fillcolor="#fce94f"];
      in_0 [label="sequence"];
      in_0 -> sra_fastq_dump;
      out_1 [label="first_read"];
      sra_fastq_dump -> out_1;
      out_2 [label="log"];
      sra_fastq_dump -> out_2;
      out_3 [label="second_read"];
      sra_fastq_dump -> out_3;
   }

**Options:**
  - **accession** (str, optional) -- Replaces accession derived from <path> in filename(s) and deflines (only for single table dump)

  - **aligned** (bool, optional) -- Dump only aligned sequences

  - **aligned-region** (str, optional) -- Filter by position on genome. Name can either be accession.version (ex:NC_000001.10) or file specific name (ex:"chr1" or "1"). "from" and "to" are 1-based coordinates. <name[:from-to]>

  - **clip** (bool, optional) -- Apply left and right clips

  - **dd-blocksize** (str, optional)    - default value: 256k

  - **defline-qual** (str, optional) -- Defline format specification for quality.

  - **defline-seq** (str, optional) -- Defline format specification for sequence.

  - **disable-multithreading** (bool, optional) -- disable multithreading

  - **dumpbase** (bool, optional) -- Formats sequence using base space (default for other than SOLiD).

  - **dumpcs** (bool, optional) -- Formats sequence using color space (default for SOLiD),"cskey" may be specified for translation.

  - **fasta** (int, optional) -- FASTA only, no qualities, optional line wrap width (set to zero for no wrapping). <[line width]>

  - **helicos** (bool, optional) -- Helicos style defline

  - **legacy-report** (bool, optional) -- use legacy style "Written spots" for tool

  - **log-level** (str, optional) -- Logging level as number or enum string One of (fatal|sys|int|err|warn|info) or (0-5). Current/default is warn. <level>

  - **matepair-distance** (str, optional) -- Filter by distance beiween matepairs. Use "unknown" to find matepairs split between the references. Use from-to to limit matepair distance on the same reference. <from-to|unknown>

  - **maxSpotId** (int, optional) -- Maximum spot id to be dumped. Use with "minSpotId" to dump a range.

  - **max_cores** (int, optional) -- Maximum number of cores available on the cluster
    - default value: 10

  - **minReadLen** (int, optional) -- Filter by sequence length >= <len>

  - **minSpotId** (int, optional) -- Minimum spot id to be dumped. Use with "maxSpotId" to dump a range.

  - **ncbi_error_report** (str, optional) -- Control program execution environment report generation (if implemented). One of (never|error|always). Default is error. <error>

  - **offset** (int, optional) -- Offset to use for quality conversion, default is 33

  - **origfmt** (bool, optional) -- Defline contains only original sequence name

  - **qual-filter** (bool, optional) -- Filter used in early 1000 Genomes data: no sequences starting or ending with >= 10N

  - **qual-filter-1** (bool, optional) -- Filter used in current 1000 Genomes data

  - **read-filter** (str, optional) -- Split into files by READ_FILTER value optionally filter by value: pass|reject|criteria|redacted

  - **readids** (bool, optional) -- Append read id after spot id as "accession.spot.readid" on defline.

  - **skip-technical** (bool, optional) -- Dump only biological reads

  - **split-spot** (str, optional) -- Split spots into individual reads

  - **spot-groups** (str, optional) -- Filter by SPOT_GROUP (member): name[,...]

  - **suppress-qual-for-cskey** (bool, optional) -- supress quality-value for cskey 

  - **table** (str, optional) -- Table name within cSRA object, default is "SEQUENCE"

  - **unaligned** (bool, optional) -- Dump only unaligned sequences

  - **verbose** (bool, optional) -- Increase the verbosity level of the program. Use multiple times for more verbosity.


**Required tools:** dd (coreutils), fastq-dump, mkfifo (coreutils), pigz

**CPU Cores:** 10

.. index:: star

star
====





**Input Connection**
  - **in/second_read**
  - **in/first_read**

**Output Connection**
  - **out/sj.out**
  - **out/log.progess**
  - **out/aligned**
  - **out/log_stdout**
  - **out/log.final**
  - **out/log_stderr**
  - **out/log.out**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      star [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> star;
      in_1 [label="second_read"];
      in_1 -> star;
      out_2 [label="aligned"];
      star -> out_2;
      out_3 [label="log.final"];
      star -> out_3;
      out_4 [label="log.out"];
      star -> out_4;
      out_5 [label="log.progess"];
      star -> out_5;
      out_6 [label="log_stderr"];
      star -> out_6;
      out_7 [label="log_stdout"];
      star -> out_7;
      out_8 [label="sj.out"];
      star -> out_8;
   }

**Options:**
  - **cores** (int, optional) -- workaround to specify cores for grid engine and threads ie
    - default value: 1

  - **genomeDir** (str, optional) -- path to the directory where genome files are stored (if runMode!=generateGenome) or will be generated (if runMode==generateGenome)

  - **readFilesCommand** (str, optional) -- command line to execute for each of the input file. This command should generate FASTA or FASTQ text and send it to stdout. For example: zcat - to uncompress .gz files, bzcat - to uncompress .bz2 files, etc.

  - **runThreadN** (int, optional) -- number of threads to run STAR
    - default value: 1


**Required tools:** rm (coreutils), star

**CPU Cores:** 1

.. index:: stringtie

stringtie
=========



    StringTie is a fast and highly efficient assembler of RNA-Seq alignments into potential
    transcripts. It uses a novel network flow algorithm as well as an optional de novo assembly step
    to assemble and quantitate full-length transcripts representing multiple splice variants for
    each gene locus. Its input can include not only the alignments of raw reads used by other
    transcript assemblers, but also alignments longer sequences that have been assembled from those
    reads.In order to identify differentially expressed genes between experiments, StringTie's
    output can be processed by specialized software like Ballgown, Cuffdiff or other programs
    (DESeq2, edgeR, etc.).

    NOTE: This step implements that part of stringtie that assembles new transcripts. If you want
    stringtie to assemble transcripts from multiple input files please use step stringtie_merge!

    https://ccb.jhu.edu/software/stringtie/


**Input Connection**
  - **in/features** (optional) Format: **['gtf', 'gff3']** - Reference assembly. Can also be passed with option G or left out for denovo assembling.
  - **in/alignments** Format: **bam**

**Output Connection**
  - **out/e2t** (optional) Format: **ctab** - Ballgown output (requires -G and -B).
  - **out/abundances** Format: **tab** - Feature abundancies (-A).
  - **out/coverage** (optional) Format: **gtf** - Coverage of the reference assmbly (-B, requires -G)
  - **out/i_data** (optional) Format: **ctab** - Ballgown output (requires -G and -B).
  - **out/features** Format: **gtf** - Contains the assempled transcripts (-o).
  - **out/e_data** (optional) Format: **ctab** - Ballgown output (requires -G and -B).
  - **out/log_stdout**
  - **out/log_stderr**
  - **out/i2t** (optional) Format: **ctab** - Ballgown output (requires -G and -B).
  - **out/t_data** (optional) Format: **ctab** - Ballgown output (requires -G and -B).


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      stringtie [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> stringtie;
      in_1 [label="features", style=filled, fillcolor="#a7a7a7"];
      in_1 -> stringtie;
      out_2 [label="abundances"];
      stringtie -> out_2;
      out_3 [label="coverage", style=filled, fillcolor="#a7a7a7"];
      stringtie -> out_3;
      out_4 [label="e2t", style=filled, fillcolor="#a7a7a7"];
      stringtie -> out_4;
      out_5 [label="e_data", style=filled, fillcolor="#a7a7a7"];
      stringtie -> out_5;
      out_6 [label="features"];
      stringtie -> out_6;
      out_7 [label="i2t", style=filled, fillcolor="#a7a7a7"];
      stringtie -> out_7;
      out_8 [label="i_data", style=filled, fillcolor="#a7a7a7"];
      stringtie -> out_8;
      out_9 [label="log_stderr"];
      stringtie -> out_9;
      out_10 [label="log_stdout"];
      stringtie -> out_10;
      out_11 [label="t_data", style=filled, fillcolor="#a7a7a7"];
      stringtie -> out_11;
   }

**Options:**
  - **B** (bool, optional) -- This switch enables the output of Ballgown input table files (\*.ctab) containing coverage data for the reference transcripts given with the -G option. (See the Ballgown documentation for a description of these files.) With this option StringTie can be used as a direct replacement of the tablemaker program included with the Ballgown distribution. The \*.ctab files will be supplied to child steps through additional connections ``out/e2t``, ``out/e_data``, ``out/i2t``, ``out/i_data`` and ``out/t_data``.

  - **G** (str, optional) -- reference annotation to use for guiding the assembly process

  - **M** (float, optional) -- Sets the maximum fraction of muliple-location-mapped reads that are allowed to be present at a given locus. Default: 0.95.

  - **a** (int, optional) -- minimum anchor length for junctions (default: 10)

  - **b** (bool, optional) -- enable output of Ballgown table files but these files will be created under the directory path given as <dir_path>

  - **c** (float, optional) -- minimum reads per bp coverage to consider for transcript assembly (default: 2.5)

  - **dd-blocksize** (str, optional) -- Provide the blocksize for dd tool.
    - default value: 2M

  - **e** (bool, optional) -- Limits the processing of read alignments to only estimate and output the assembled transcripts matching the reference transcripts given with the -G option (requires -G, recommended for -B/-b). With this option, read bundles with no reference transcripts will be entirely skipped, which may provide a considerable speed boost when the given set of reference transcripts is limited to a set of target genes, for example.

  - **f** (float, optional) -- Sets the minimum isoform abundance of the predicted transcripts as a fraction of the most abundant transcript assembled at a given locus. Lower abundance transcripts are often artifacts of incompletely spliced precursors of processed transcripts. Default: 0.1

  - **fifo** (bool, optional) -- Enable the FIFO functionality for splitting large input files.

  - **fr** (bool, optional) -- assume stranded library fr-secondstrand

  - **g** (int, optional) -- gap between read mappings triggering a new bundle (default: 50)

  - **j** (float, optional) -- minimum junction coverage (default: 1)

  - **l** (str, optional) -- Sets <label> as the prefix for the name of the output transcripts. Default: STRG

  - **m** (int, optional) -- Sets the minimum length allowed for the predicted transcripts. Default: 200

  - **p** (int, optional) -- Specify the number of processing threads (CPUs) to use for transcript assembly.
    - default value: 6

  - **rf** (bool, optional) -- assume stranded library fr-firststrand

  - **t** (bool, optional) -- disable trimming of predicted transcripts based on coverage (default: coverage trimming is enabled)

  - **v** (bool, optional) -- Turns on verbose mode, printing bundle processing details

  - **x** (str, optional) -- Ignore all read alignments (and thus do not attempt to perform transcript assembly) on the specified reference sequences. Parameter <seqid_list> can be a single reference sequence name (e.g. -x chrM) or a comma-delimited list of sequence names (e.g. -x "chrM,chrX,chrY"). This can speed up StringTie especially in the case of excluding the mitochondrial genome, whose genes may have very high coverage in some cases, even though they may be of no interest for a particular RNA-Seq analysis. The reference sequence names are case sensitive, they must match identically the names of chromosomes/contigs of the target genome against which the RNA-Seq reads were aligned in the first place.


**Required tools:** dd (coreutils), mkdir (coreutils), mkfifo (coreutils), mv (coreutils), stringtie

**CPU Cores:** 6

.. index:: stringtieMerge

stringtieMerge
==============



    # stringtie --merge <gtf.list> > outputpat/outputname

    StringTie is a fast and highly efficient assembler of RNA-Seq alignments into potential
    transcripts. merge is a mode of the StringTie tool that is used to assemble transcripts from multiple input files (assemblies). It generates a unified non-redundant set of isoforms.

    NOTE: This step implements the merging part of stringtie. If you want
    stringtie to assemble transcripts from multiple BAM files please use step stringtie!

    https://ccb.jhu.edu/software/stringtie/

**Input Connection**
  - **in/reference** (optional) Format: **['gtf', 'gff3']** - Reference assembly. Can also be passed with option G or left out for denovo assembling.
  - **in/features** Format: **['gtf', 'gff3']** - Feature annotations to be merged.

**Output Connection**
  - **out/assemblies**
  - **out/features** Format: **gtf**
  - **out/log_stderr**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      stringtieMerge [style=filled, fillcolor="#fce94f"];
      in_0 [label="features"];
      in_0 -> stringtieMerge;
      in_1 [label="reference", style=filled, fillcolor="#a7a7a7"];
      in_1 -> stringtieMerge;
      out_2 [label="assemblies"];
      stringtieMerge -> out_2;
      out_3 [label="features"];
      stringtieMerge -> out_3;
      out_4 [label="log_stderr"];
      stringtieMerge -> out_4;
   }

**Options:**
  - **F** (float, optional) -- minimum input transcript FPKM to include in the merge (default: 1.0)

  - **G** (str, optional) -- reference annotation to include in the merging (GTF/GFF3)

  - **T** (float, optional) -- minimum input transcript TPM to include in the merge (default: 1.0)

  - **c** (int, optional) -- minimum input transcript coverage to include in the merge (default: 0)

  - **f** (float, optional) -- minimum isoform fraction (default: 0.01)

  - **g** (int, optional) -- gap between transcripts to merge together (default: 250)

  - **i** (bool, optional) -- keep merged transcripts with retained introns; by default

  - **l** (str, optional) -- name prefix for output transcripts (default: MSTRG)

  - **m** (int, optional) -- minimum input transcript length to include in the merge (default: 50)

  - **output_prefix** (str, optional) -- Prefix used in the utput directory.
    - default value: merge

  - **p** (int, optional) -- Number of cores
    - default value: 2


**Required tools:** mkdir (coreutils), mv (coreutils), printf (coreutils), stringtie

**CPU Cores:** 2

.. index:: stringtie_prepDE

stringtie_prepDE
================


StringTie is a fast and highly efficient assembler of RNA-Seq alignments into potential
    transcripts. prepDE.py is Python script to extract this read count information directly from the
    files generated by StringTie (run with the -e parameter). It generates two CSV files containing
    the count matrices for genes and transcripts, using the coverage values found in the output of
    stringtie -e

    Here, all the transcripts.gtf files from a previous stringtie call are collected, written to a
    file and provided to the prepDE.py script for conversion (all at once).

    NOTE: This step implements the prepDE.py part of stringtie. If you want stringtie to assemble
    transcripts from multiple BAM files please or merge assemblies use step stringtie or
    stringtie_merge, resp.!

    https://ccb.jhu.edu/software/stringtie/


**Input Connection**
  - **in/features**

**Output Connection**
  - **out/legend**
  - **out/log_stdout**
  - **out/log_stderr**
  - **out/transcript_matrix**
  - **out/gene_matrix**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      stringtie_prepDE [style=filled, fillcolor="#fce94f"];
      in_0 [label="features"];
      in_0 -> stringtie_prepDE;
      out_1 [label="gene_matrix"];
      stringtie_prepDE -> out_1;
      out_2 [label="legend"];
      stringtie_prepDE -> out_2;
      out_3 [label="log_stderr"];
      stringtie_prepDE -> out_3;
      out_4 [label="log_stdout"];
      stringtie_prepDE -> out_4;
      out_5 [label="transcript_matrix"];
      stringtie_prepDE -> out_5;
   }

**Options:**
  - **cluster** (bool, optional) -- whether to cluster genes that overlap with different gene IDs, ignoring ones with geneID pattern (see below)

  - **key** (str, optional) -- if clustering, what prefix to use for geneIDs assigned by this script [default: prepG]

  - **length** (int, optional) -- the average read length [default: 75]

  - **pattern** (str, optional) -- a regular expression that selects the sample subdirectories

  - **run_id** (str, optional) -- A name for the run. Since this step merges multiple samples into a single one, the run_id cannot be the sample name anymore.
    - default value: prepDEall

  - **string** (str, optional) -- if a different prefix is used for geneIDs assigned by StringTie [default: MSTRG


**Required tools:** prepDE, printf (coreutils)

**CPU Cores:** 6

.. index:: subsetMappedReads

subsetMappedReads
=================



    subsetMappedReads selects a provided number of mapped reads from a
    file in .sam or .bam format. Depending on the set options the
    first N mapped reads and their mates (for paired end sequencing)
    are returned in .sam format. If the number of requested reads
    exceeds the number of available mapped reads, all mapped reads are
    returned.

**Input Connection**
  - **in/alignments**

**Output Connection**
  - **out/alignments**
  - **out/log**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      subsetMappedReads [style=filled, fillcolor="#fce94f"];
      in_0 [label="alignments"];
      in_0 -> subsetMappedReads;
      out_1 [label="alignments"];
      subsetMappedReads -> out_1;
      out_2 [label="log"];
      subsetMappedReads -> out_2;
   }

**Options:**
  - **Nreads** (str, required) -- Number of reads to extract from input file. 

  - **genome-faidx** (str, required)
  - **paired_end** (bool, required) -- The reads are expected to have a mate, due to paired end sequencing.


**Required tools:** cat (coreutils), dd (coreutils), head (coreutils), pigz, samtools

**CPU Cores:** 1

.. index:: tcount2gcount

tcount2gcount
=============



    This step converts transcript based count files, e.g., from kallosto or
    salmon, into a gene base count files of the same format by summing
    counts of transcripts of the same gene.

**Input Connection**
  - **in/annotation**
  - **in/counts**

**Output Connection**
  - **out/counts**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      tcount2gcount [style=filled, fillcolor="#fce94f"];
      in_0 [label="annotation"];
      in_0 -> tcount2gcount;
      in_1 [label="counts"];
      in_1 -> tcount2gcount;
      out_2 [label="counts"];
      tcount2gcount -> out_2;
   }

**Options:**
  - **cores** (int, optional) -- workaround to specify cores for grid engine and threads ie
    - default value: 1

  - **kallisto-extended** (bool, optional) -- writes extended format includign tpm. 

  - **m** (str, optional) -- transcript to gene mapping file. Required Format example (per row): ENST00000527779.1 ENSG00000137692.11 or gtf

  - **t** (str, optional) -- source tool of input file. Possible values: kallisto, salmon


**CPU Cores:** 1

.. index:: tophat2

tophat2
=======



    TopHat is a fast splice junction mapper for RNA-Seq reads.
    It aligns RNA-Seq reads to mammalian-sized genomes using the ultra
    high-throughput short read aligner Bowtie, and then analyzes the mapping
    results to identify splice junctions between exons.

    http://tophat.cbcb.umd.edu/

    typical command line::

        tophat [options]* <index_base> <reads1_1[,...,readsN_1]>         [reads1_2,...readsN_2]

    Tested on release: TopHat v2.0.13

**Input Connection**
  - **in/second_read**
  - **in/first_read**

**Output Connection**
  - **out/prep_reads**
  - **out/deletions**
  - **out/unmapped**
  - **out/misc_logs**
  - **out/alignments**
  - **out/log_stderr**
  - **out/align_summary**
  - **out/insertions**
  - **out/junctions**


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
  - **bowtie1** (bool, optional) -- Use bowtie1. Default: bowtie2

  - **color** (bool, optional) -- Solid - color space

  - **color-out** (bool, optional) -- Colored output

  - **index** (str, required) -- Path to genome index for tophat2

  - **integer-quals** (bool, optional) -- Provide/Use (?) integer qualities.

  - **library_type** (str, required) -- The default is unstranded (fr-unstranded). If either fr-firststrand or fr-secondstrand is specified, every read alignment will have an XS attribute tag as explained below. Consider supplying library type options below to select the correct RNA-seq protocol.(https://ccb.jhu.edu/software/tophat/manual.shtml)
    - possible values: 'fr-unstranded', 'fr-firststrand', 'fr-secondstrand'


  - **max-deletion-length** (int, optional) -- Max size of deletion

  - **max-insertion-length** (int, optional) -- Max size of insertion

  - **max-intron-length** (int, optional) -- maximal intron length

  - **max-multihits** (int, optional) -- Maximal number of multiple hits

  - **min-anchor** (int, optional) -- Size of minimal anchor.

  - **min-intron-length** (int, optional) -- Minimal intron length

  - **phred64-quals** (bool, optional) -- Qualities are phred64 qualities (same as solexa1.3-quals).

  - **prefilter-multihits** (bool, optional) -- for -G/--GTF option, enable an initial bowtie search against the genome

  - **quals** (bool, optional) -- Provide/Use (?) qualities.

  - **read-edit-dist** (int, optional) -- Read edit distance

  - **read-gap-length** (int, optional) -- Size of gap length

  - **read-mismatches** (int, optional)
  - **read-realign-edit-dist** (int, optional) -- Read alignment distance. Default: read-edit-dist + 1.

  - **solexa-quals** (bool, optional) -- Qualities are solexa qualities.

  - **solexa1.3-quals** (bool, optional) -- Qualities are solexa1.3 qualities (same as phred64-quals).

  - **splice-mismatches** (int, optional) -- Number of splice mismatches
    - possible values: '0', '1', '2'


  - **supress-hits** (bool, optional) -- Supress hits

  - **transcriptome-max-hits** (int, optional) -- Max hits in transcriptome


**Required tools:** mkdir (coreutils), mv (coreutils), tar, tophat2

**CPU Cores:** 6

.. index:: trim_galore

trim_galore
===========



    A wrapper tool around Cutadapt and FastQC to consistently apply quality and adapter trimming to
    FastQ files, with some extra functionality for MspI-digested RRBS-type (Reduced Representation
    Bisufite-Seq) libraries.

    https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/

    Note for RRBS using the NuGEN Ovation RRBS System 1-16 kit:

    Owing to the fact that the NuGEN Ovation kit attaches a varying number of nucleotides (0-3) after each MspI
    site Trim Galore should be run WITHOUT the option --rrbs. This trimming is accomplished in a subsequent
    diversity trimming step afterwards (see their manual).

    Note for RRBS using MseI:

    If your DNA material was digested with MseI (recognition motif: TTAA) instead of MspI it is NOT necessary
    to specify --rrbs or --non_directional since virtually all reads should start with the sequence
    'TAA', and this holds true for both directional and non-directional libraries. As the end-repair of 'TAA'
    restricted sites does not involve any cytosines it does not need to be treated especially. Instead, simply
    run Trim Galore! in the standard (i.e. non-RRBS) mode.


**Input Connection**
  - **in/second_read**
  - **in/first_read**

**Output Connection**
  - **out/second_read**
  - **out/first_read_report**
  - **out/stderr**
  - **out/second_read_fastqc_zip**
  - **out/first_read**
  - **out/stdout**
  - **out/second_read_fastqc_html**
  - **out/first_read_fastqc_html**
  - **out/first_read_fastqc_zip**
  - **out/second_read_report**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      trim_galore [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> trim_galore;
      in_1 [label="second_read"];
      in_1 -> trim_galore;
      out_2 [label="first_read"];
      trim_galore -> out_2;
      out_3 [label="first_read_fastqc_html"];
      trim_galore -> out_3;
      out_4 [label="first_read_fastqc_zip"];
      trim_galore -> out_4;
      out_5 [label="first_read_report"];
      trim_galore -> out_5;
      out_6 [label="second_read"];
      trim_galore -> out_6;
      out_7 [label="second_read_fastqc_html"];
      trim_galore -> out_7;
      out_8 [label="second_read_fastqc_zip"];
      trim_galore -> out_8;
      out_9 [label="second_read_report"];
      trim_galore -> out_9;
      out_10 [label="stderr"];
      trim_galore -> out_10;
      out_11 [label="stdout"];
      trim_galore -> out_11;
   }

**Options:**
  - **adapter** (str, optional) -- Adapter sequence to be trimmed. If not specified explicitly, Trim Galore willtry to auto-detectwhether the Illumina universal, Nextera transposase or Illuminasmall RNA adapter sequence was used.Also see '--illumina', '--nextera' and '--small_rna'. If no adapter can be detected within thefirst 1 million sequencesof the first file specified Trim Galore defaults to '--illumina'.

  - **adapter2** (str, optional) -- Optional adapter sequence to be trimmed off read 2 of paired-end files. Thisoption requires'--paired' to be specified as well. If the libraries to be trimmedare smallRNA then a2 will be setto the Illumina small RNA 5' adapter automatically(GATCGTCGGACT).

  - **clip_R1** (int, optional) -- Instructs Trim Galore to remove <int> bp from the 5' end of read 1 (or single-endreads). This maybe useful if the qualities were very poor, or if there is somesort of unwanted bias at the 5' end.Default: OFF.

  - **clip_R2** (int, optional) -- Instructs Trim Galore to remove <int> bp from the 5' end of read 2 (paired-end readsonly). Thismay be useful if the qualities were very poor, or if there is some sortof unwanted bias at the 5'end. For paired-end BS-Seq, it is recommended to removethe first few bp because the end-repairreaction may introduce a bias towards lowmethylation. Please refer to the M-bias plot section in theBismark User Guide forsome examples. Default: OFF.

  - **dd-blocksize** (str, optional)    - default value: 2M

  - **dont_gzip** (bool, optional) -- Output files won't be compressed with GZIP. This option overrides --gzip.

  - **e** (float, optional) -- Maximum allowed error rate (no. of errors divided by the length of the matchingregion) (default:0.1)

  - **fastqc** (bool, optional) -- Run FastQC in the default mode on the FastQ file once trimming is complete.

  - **fastqc_args** (str, optional) -- Passes extra arguments to FastQC. If more than one argument is to be passedto FastQC they must bein the form "arg1 arg2 etc.". An example would be:--fastqc_args "--nogroup --outdir /home/". Passingextra arguments willautomatically invoke FastQC, so --fastqc does not have to bespecifiedseparately.

  - **first_read** (str, optional) -- Part of the file name that marks all files containing sequencing data of the first read. Example: '_R1' or '_1'. Default: 'R1'
    - default value: R1

  - **gzip** (bool, optional) -- Compress the output file with GZIP. If the input files are GZIP-compressedthe output files willautomatically be GZIP compressed as well. As of v0.2.8 the compression will take place on the fly.

  - **illumina** (bool, optional) -- Adapter sequence to be trimmed is the first 13bp of the Illumina universal adapter'AGATCGGAAGAGC' instead of the default auto-detection of adapter sequence.

  - **keep** (bool, optional) -- Keep the quality trimmed intermediate file. Default: off, which meansthe temporary file is beingdeleted after adapter trimming. Only hasan effect for RRBS samples since other FastQ files are nottrimmedfor poor qualities separately.

  - **length** (int, optional) -- Discard reads that became shorter than length INT because of eitherquality or adapter trimming. Avalue of '0' effectively disablesthis behaviour. Default: 20 bp.For paired-end files, both readsof a read-pair need to be longer than<INT> bp to be printed out to validated paired-end files (seeoption --paired).If only one read became too short there is the possibility of keeping suchunpairedsingle-end reads (see --retain_unpaired). Default pair-cutoff: 20 bp.

  - **length_1** (int, optional) -- Unpaired single-end read length cutoff needed for read 1 to be written to '.unpaired_1.fq 'output file. These reads may be mapped in single-end mode.Default: 35 bp.

  - **length_2** (int, optional) -- Unpaired single-end read length cutoff needed for read 2 to be written to '.unpaired_2.fq 'output file. These reads may be mapped in single-end mode.Default: 35 bp.

  - **max_length** (int, optional) -- Discard reads that are longer than <INT> bp after trimming. This is only advised for smallRNAsequencing to remove non-small RNA sequences.

  - **max_n** (int, optional) -- The total number of Ns (as integer) a read may contain before it will be removed altogether.In apaired-end setting, either read exceeding this limit will result in the entirepair being removedfrom the trimmed output files.

  - **nextera** (bool, optional) -- Adapter sequence to be trimmed is the first 12bp of the Nextera adapter 'CTGTCTCTTATA' instead ofthe default auto-detection of adapter sequence.

  - **non_directional** (bool, optional) -- Selecting this option for non-directional RRBS libraries will screenquality-trimmed sequences for'CAA' or 'CGA' at the start of the readand, if found, removes the first two basepairs. Like withthe option '--rrbs' this avoids using cytosine positions that were filled-induring the end-repairstep. '--non_directional' requires '--rrbs' tobe specified as well. Note that this option doesnot set '--clip_r2 2' inpaired-end mode.

  - **paired** (bool, optional) -- This option performs length trimming of quality/adapter/RRBS trimmed reads forpaired-end files. Topass the validation test, both sequences of a sequence pairare required to have a certain minimumlength which is governed by the option--length (see above). If only one read passes this lengththreshold theother read can be rescued (see option --retain_unpaired). Using this option letsyoudiscard too short read pairs without disturbing the sequence-by-sequence orderof FastQ files whichis required by many aligners.Trim Galore! expects paired-end files to be supplied in a pairwisefashion, e.g.file1_1.fq file1_2.fq SRR2_1.fq.gz SRR2_2.fq.gz ... .

  - **phred33** (bool, optional) -- Instructs Cutadapt to use ASCII+33 quality scores as Phred scores (Sanger/Illumina 1.9+ encoding)for quality trimming. Default: ON.

  - **phred64** (bool, optional) -- Instructs Cutadapt to use ASCII+64 quality scores as Phred scores(Illumina 1.5 encoding) forquality trimming.

  - **pigz-blocksize** (str, optional)    - default value: 2048

  - **quality** (int, optional) -- Trim low-quality ends from reads in addition to adapter removal. For RRBS samples, quality trimmingwill be performed first, and adaptertrimming is carried in a second round. Other files are qualityand adaptertrimmed in a single pass. The algorithm is the same as the one used by BWA(Subtract INTfrom all qualities; compute partial sums from all indicesto the end of the sequence; cut sequence atthe index at which the sum isminimal). Default Phred score: 20.

  - **retain_unpaired** (bool, optional) -- If only one of the two paired-end reads became too short, the longerread will be written to either'.unpaired_1.fq' or '.unpaired_2.fq' output files. The length cutoff for unpaired single-endreads isgoverned by the parameters -r1/--length_1 and -r2/--length_2. Default: OFF.

  - **rrbs** (bool, optional) -- Specifies that the input file was an MspI digested RRBS sample (recognitionsite: CCGG). Single-endor Read 1 sequences (paired-end) which were adapter-trimmedwill have a further 2 bp removed fromtheir 3' end. Sequences which were merelytrimmed because of poor quality will not be shortenedfurther. Read 2 of paired-endlibraries will in addition have the first 2 bp removed from the 5' end(by setting '--clip_r2 2'). This is to avoid using artificial methylation calls from thefilled-incytosine positions close to the 3' MspI site in sequenced fragments. This option is notrecommended for users of the NuGEN ovation RRBS System 1-16kit (see below).

  - **second_read** (str, optional) -- Part of the file name that marks all files containing sequencing data of the second read. Example: '_R2' or '_2'. Default: 'R2'
    - default value: R2

  - **small_rna** (bool, optional) -- Adapter sequence to be trimmed is the first 12bp of the Illumina Small RNA 3' Adapter'TGGAATTCTCGG' instead of the default auto-detection of adapter sequence. Selectingto trimsmallRNA adapters will also lower the --length value to 18bp. If the smallRNAlibraries arepaired-end then a2 will be set to the Illumina small RNA 5' adapterautomatically (GATCGTCGGACT)unless -a 2 had been defined explicitly.

  - **stringency** (int, optional) -- Overlap with adapter sequence required to trim a sequence. Defaults to avery stringent setting of1, i.e. even a single bp of overlapping sequencewill be trimmed off from the 3' end of any read.

  - **three_prime_clip_R1** (int, optional) -- Instructs Trim Galore to remove <int> bp from the 3' end of read 1 (or single-endreads) AFTERadapter/quality trimming has been performed. This may remove some unwantedbias from the 3' end thatis not directly related to adapter sequence or basecall quality.Default: OFF.

  - **three_prime_clip_R2** (int, optional) -- Instructs Trim Galore to remove <int> bp from the 3' end of read 2 AFTERadapter/quality trimminghas been performed. This may remove some unwanted bias fromthe 3' end that is not directly relatedto adapter sequence or basecall quality.Default: OFF.

  - **trim-n** (bool, optional) -- Removes Ns from either side of the read. This option does currently not work in RRBS mode.

  - **trim1** (bool, optional) -- Trims 1 bp off every read from its 3' end. This may be needed for FastQ files thatare to bealigned as paired-end data with Bowtie. This is because Bowtie (1) regardsalignments like this: R1---------------------------> or this: -----------------------> R1 R2 <---------------------------<----------------- R2as invalid (whenever a start/end coordinate is contained within the otherread).NOTE: If you are planning to use Bowtie2, BWA etc. you don't need to specify this option.


**Required tools:** cutadapt, trim_galore

**CPU Cores:** 4

.. index:: trimmomatic

trimmomatic
===========


**Input Connection**
  - **in/second_read**
  - **in/first_read**

**Output Connection**
  - **out/forward**
  - **out/reverse**
  - **out/reverse.unpaired**
  - **out/log_stderr**
  - **out/forward.unpaired**
  - **out/log**


.. graphviz::

   digraph foo {
      rankdir = LR;
      splines = true;
      graph [fontname = Helvetica, fontsize = 12, size = "14, 11", nodesep = 0.2, ranksep = 0.3];
      node [fontname = Helvetica, fontsize = 12, shape = rect];
      edge [fontname = Helvetica, fontsize = 12];
      trimmomatic [style=filled, fillcolor="#fce94f"];
      in_0 [label="first_read"];
      in_0 -> trimmomatic;
      in_1 [label="second_read"];
      in_1 -> trimmomatic;
      out_2 [label="forward"];
      trimmomatic -> out_2;
      out_3 [label="forward.unpaired"];
      trimmomatic -> out_3;
      out_4 [label="log"];
      trimmomatic -> out_4;
      out_5 [label="log_stderr"];
      trimmomatic -> out_5;
      out_6 [label="reverse"];
      trimmomatic -> out_6;
      out_7 [label="reverse.unpaired"];
      trimmomatic -> out_7;
   }

**Options:**
  - **base_file_name** (int, optional) -- The prefix commoon to all output files, replacing the run-ID

  - **jar_path** (str, optional) -- Path to trimmomatic.jar
    - default value: trimmomatic

  - **phred-base** (int, optional) -- PHRED-base
    - possible values: '33', '64'


  - **steps** (list, required) -- List defining the analysis, in terms of trimmomatic-commands

  - **threads** (int, optional) -- Number of threads to use
    - default value: 1


**Required tools:** java

**CPU Cores:** 1

