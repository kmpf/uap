..
  This is the documentation for uap. Please keep lines under
  80 characters if you can and start each sentence on a new line as it
  decreases maintenance and makes diffs more readable.

.. title:: Configuration of uap

..
  This document aims to describe how to configure **uap**.

.. _analysis_configuration:

***************************
Analysis Configuration File
***************************

**uap** requires a |yaml_link| file which contains all information
about the data analysis.
These files are called configuration files.

A configuration file describes a complete analysis.
Configurations consist of four sections (let's just call them sections,
although technically, they are keys):

**Mandatory Sections**

  * ``destination_path`` -- points to the directory where the result files,
    annotations and temporary files are written to
  * ``constants`` -- defines constants for later use (define repeatedly used
    values as constants to increase readability of the following sections)
  * ``steps`` -- defines the source and processing steps and their order
  * ``lmod`` -- if lmod is used paths can be specified here to ignor user env
  * ``tools`` -- defines all tools used in the analysis and how to determine
    their versions (for later reference)

**Optional Sections**

  * ``cluster`` -- if **uap** is required to run on a HPC cluster some default
    parameters can be set her

Please refer to the |yaml_link| definition for the correct notation used in
that file.

Sections of a Configuration File
================================

.. _config-file-destination-path:

``destination_path`` Section
----------------------------

The value of ``destination_path`` is the directory where **uap** is going
to store the created files.

.. It is possible to use a different directory for volatile files (see ).

.. code-block:: yaml

    destination_path: "/path/to/workflow/output"


``base_working_directory`` Section
----------------------------------

The value of ``base_working_directory`` is the directory where **uap**
changes to before declaring all steps and it defaults to the location
of the configuration file. All configured paths can be set relatively
to this directory.

.. code-block:: yaml

    base_working_directory: "/path/to/workflow/output"


``constants`` Section
---------------------

This section is the place where repeatedly used constants should be defined.
For instance absolute paths to the genome index files can be defined as
constant.

.. code-block:: yaml

   - &genome_faidx
        genomes/animalia/chordata/mammalia/primates/homo_sapiens/hg19/samtools_faidx/hg19_all_chr_UCSC-download-B7ceRp9K/hg19_all_chr_UCSC-download.fasta.fai

Later on the value can be reused by typing ``*genome_faidx``.
**There are no restrictions about what can be defined here.**

.. _config-file-steps:

``steps`` Section
-----------------

The ``steps`` section is the core of the analysis file, because it defines when
steps are executed and how they depend on each other.
All available steps are described in detail in the steps documentation:
:doc:`steps`.
The ``steps`` section contains an entry (technically a key) for every step.
Every step name **must** be unique.

.. note::

   Please be aware that the |pyyaml_link|, the YAML parser used by uap, does not
   complain about keys with the same name.
   But drops one of the duplicates without giving an error.

There are two ways to name a step to allow multiple steps of the same type and
still ensure unique naming:

.. code-block:: yaml

    steps:
        # here, the step name is unchanged, it's a cutadapt step which is also
        # called 'cutadapt'
        cutadapt:
            ... # options following

        # here, we also insert a cutadapt step, but we give it a different name:
        # 'clip_adapters'
        clip_adapters (cutadapt):
            ... # options following


Now let's have a look at the two different types of steps which constitute
an **uap** analaysis.

.. _config_file_source_steps:

Source Steps
^^^^^^^^^^^^

Source steps are the only steps which are allowed to use or create data
outside the ``destination_path``.
Feature of source steps:

* they provide the input files for the following steps
* they can start processes e.g. to download files or demultiplex reads
* they do not depend on previous steps
* they are the root nodes of the analysis graph

If you want to work with fastq files, you should use the ``fastq_source``
step to import the required files.
Such a step definition would look like this:

.. code-block:: yaml

    steps:
        input_step (fastq_source):
        pattern: /Path/to/fastq-files/*.gz
        group: ([SL]\w+)_R[12]-00[12].fastq.gz
        sample_id_prefix: MyPrefix
        first_read: '_R1'
        second_read: '_R2'
        paired_end: True

The options of the ``fastq_source`` step are described at :doc:`steps`.
The ``group`` option takes a regular expression (regexp).
You can test your regular expression at |pythex_link|.

.. _config_file_processing_steps:

Processing Steps
^^^^^^^^^^^^^^^^

Processing steps depend upon one or more preceding steps.
They use their output files and process them.
Output files of processing steps are automatically named and saved by **uap**.
A complete list of available options per step can be found at :doc:`steps`
or by using the :ref:`uap-steps`.

.. _config_file_keywords:

Reserved Keywords for Steps
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _config_file_depends:

**_depends:**

  Dependencies are defined via the ``_depends`` key which may either be ``null``,
  a step name, or a list of step names.

.. code-block:: yaml

    steps:
        # the source step which depends on nothing
        fastq_source:
            # ...

        run_folder_source:
            # ...

        # the first processing step, which depends on the source step
        cutadapt:
            _depends: [fastq_source, run_folder_source]

        # the second processing step, which depends on the cutadapt step
        fix_cutadapt:
            _depends: cutadapt

.. _config_file_connect:

**_connect:**

  Normally steps connected with ``_depends`` do pass data along by defining
  so called connections.
  If the name of an output connection matches the name of an input connection
  of a succeeding step the data gets passed on automatically.
  But, sometimes the user wants to force the connection of differently named
  connections.
  This can be done with the ``_connect`` keyword.
  A common usage is to connect downloaded data with a
  :ref:`config_file_processing_steps`.

.. code-block:: yaml

    steps:
        # Source step to download i.e. sequence of chr1 of some species
        chr1 (raw_url_source):
            ...

        # Download chr2 sequence
        chr2 (raw_url_source):
            ...

        merge_fasta_files:
            _connect:
                in/sequence:
                    - chr1/raw
                    - chr2/raw
            # Equivalent to:
            # _connect:
            #     in/sequence: [chr1/raw, chr2/raw]

  The examples shows how the ``raw_url_source`` output connection ``raw`` is
  connected to the input connection ``sequence`` of the ``merge_fasta_files``
  step.

.. _config_file_break:

**_BREAK:**

  If you want to cut off entire branches of the step graph, set the ``_BREAK``
  flag in a step definition, which will force the step to produce no runs
  (which will in turn give all following steps nothing to do, thereby
  effectively disabling these steps):


.. code-block:: yaml

    steps:
        fastq_source:
            # ...

        cutadapt:
            _depends: fastq_source

        # this step and all following steps will not be executed
        fix_cutadapt:
            _depends: cutadapt
            _BREAK: true

.. _config_file_volatile:

**_volatile:**

  Steps can be marked with ``_volatile: yes``.
  This flag tells **uap** that the output files of the marked step are only
  intermediate results.

.. code-block:: yaml

    steps:
        # the source step which depends on nothing
        fastq_source:
            # ...

        # this steps output can be deleted if all depending steps are finished
        cutadapt:
            _depends: fastq_source
            _volatile: yes
            # same as:
            # _volatile: True

        # if fix_cutadapt is finished the output files of cutadapt can be
        # volatilized
        fix_cutadapt:
            _depends: cutadapt

If all steps depending on the intermediate step are finished **uap** tells the
user that he can free disk space.
The message is output if the :ref:`status <uap-status>` is checked and looks
like this::

   Hint: You could save 156.9 GB of disk space by volatilizing 104 output files.
   Call 'uap <project-config>.yaml volatilize --srsly' to purge the files.

**uap** is going to replace the output files by placeholder files if the user
executes the :ref:`volatilize <uap-volatilize>` command.

.. _config_file_cluster_submit_options:

**_cluster_submit_options**

    This string contains the entire submit options which will be set in the
    submit script.
    This option allows to overwrite the values set in
    :ref:`default_submit_options <config_file_default_submit_options>`.

.. _config_file_cluster_pre_job_command:

**_cluster_pre_job_command**

    This string contains command(s) that are executed **BEFORE uap** is started
    on the cluster.
    This option allows to overwrite the values set in
    :ref:`default_pre_job_command <config_file_default_pre_job_command>`.

.. _config_file_cluster_post_job_command:

**_cluster_post_job_command**

    This string contains command(s) that are executed **AFTER uap** did finish
    on the cluster.
    This option allows to overwrite the values set in
    :ref:`default_post_job_command <config_file_default_post_job_command>`.

.. _config_file_cluster_job_quota:

**_cluster_job_quota**

    This option defines the number of jobs of the same type that can
    run simultaneously on a cluster.
    This option allows to overwrite the values set in
    :ref:`default_job_quota <config_file_default_job_quota>`.

.. _config_file_tools:

### `tools` Section

The ``tools`` section lists all programs required for the execution of a
particular analysis.
An example tool configuration looks like this:

.. code-block:: yaml

   tools:

        # you don't have to specify a path if the tool can be found in $PATH
        cat:
            path: cat
            get_version: --version

        # you have to specify a path if the tool can not be found in $PATH
        some-tool:
            path: /path/to/some-tool
            get_version: --version

       # if the output is not sesetive to the tool version it can be ignored
       mv:
          ignore_version: True

       pigz:
           path: pigz
           get_version: --version
           exit_code: 0


**uap** uses the ``path``, ``get_version``, and ``exit_code`` information to
control the availability of a tool.
This is particularly useful on cluster systems were software can be dynamically
loaded and unloaded.
**uap** logs the version of every used tool.
If ``get_version`` and ``exit_code`` is not set, **uap** tries to determine the
version by calling the program without command-line arguments.
``get_version`` is the command line argument (e.g. ``--version``) required to
get the version information.
``exit_code`` is the value returned by ``echo $?`` after trying to determine
the version e.g. by running ``pigz --version``.
If not set ``exit_code`` defaults to 0, ``get_version`` to ``--version``,
``ignore_version`` to ``False`` and ``path`` to the tool name.

Some tools are configured by default. Theire configuration will be logged
in the result annotation but they do not have to be made explicitly in the
configuration yaml. These are tools that come with the UAP installation
in ``<UAP path>/tools`` and these |coreutils|: basename, cat, cp, cut, date,
dd, dirname, du, head, ln, ls, mkdir, mkfifo, mv, paste, printf, pwd, seq,
sleep, sort, rm, tail, tee, tr, uniq, wc. The ``ignore_version`` of these
|coreutils| defaults to ``True``.

To use |lmod_link| to load an unload a tool you can specify the
``module_name`` option:

.. code-block:: yaml

   tools:

       pigz:
           path: pigz
           get_version: --version
           exit_code: 0
           module_name: pigz/version


.. _config_file_lmod:

``lmod`` Section
-------------------

This section is optional and specifies the |lmod_link| utility. It is
only required if |lmod_link| is not loaded and ``module_name`` is
used in the ``tools`` section.

.. code-block:: yaml

    lmod:
        path: /path/to/lmod/executable
        module_path: /colon/seperated/paths/to/the/used/modules

``path`` defaults to ``$LMOD_CMD`` and ``module_path`` to ``$MODULEPATH``
of the user environment.


.. _config_file_cluster:

### ``cluster`` Section
-------------------

The ``cluster`` section is required only if the analysis is executed on a
system using a cluster engine like |uge_link| or |slurm_link|.
This section interacts tightly with the
An example ``cluster`` section looks like this:

.. code-block:: yaml

    cluster:
        default_submit_options: "-pe smp #{CORES} -cwd -S /bin/bash -m as -M me@example.com -l h_rt=1:00:00 -l h_vmem=2G"
        default_pre_job_command: "echo 'Started the run!'"
        default_post_job_command: "echo 'Finished the run!'"
        default_job_quota: 5

.. _config_file_default_submit_options:

**default_submit_options**

    This is the default submit options string which replaces the
    :ref:`#{SUBMIT_OPTIONS} <submit_template_submit_options>` placeholder in
    the :ref:`submit script template <submit_template>`.
    It is **mandatory** to set this value.

.. _config_file_default_pre_job_command:

**default_pre_job_command**

    This string contains the default commands which will be executed
    **BEFORE uap** is started on the cluster.
    It will replace the
    :ref:`#{PRE_JOB_COMMAND} <submit_template_pre_job_command>` placeholder in
    the :ref:`submit script template <submit_template>`.
    If mutliple commands shall be executed separate those with ``\n``.
    It is **optional** to set this value.

.. _config_file_default_post_job_command:

**default_post_job_command**

    This string contains the default commands which will be executed
    **AFTER uap** is started on the cluster.
    It will replace the
    :ref:`#{POST_JOB_COMMAND} <submit_template_post_job_command>` placeholder in
    the :ref:`submit script template <submit_template>`.
    If mutliple commands shall be executed separate those with ``\n``.
    It is **optional** to set this value.

.. _config_file_default_job_quota:

**default_job_quota:**

    This option defines the number of jobs of the same type that can
    run simultaneously on a cluster.
    A value *0* means no limit is applied.
    It is **optional** to set this value, if the value is not provided it
    defaults to *0*.

Example Configurations
======================

Example configurations can be found in **uap**'s ``example-configurations``
folder.
More information about these examples can be found in :doc:`how-to`.

.. _cluster_configuration:

**************************
# Cluster Configuration File
**************************

The cluster configuration file resides at::

    $ ls -la $(dirname $(which uap))/cluster/cluster-specific-commands.yaml

This YAML file contains a dictionary for every cluster type.
An example file is shown here:

.. code-block:: yaml

   # Configuration for a UGE cluster engine
   uge:
       # Command to get version information
       identity_test: ['qstat', '-help']
       # The expected output of identity_test for this cluster engine
       identity_answer: 'UGE'
       # Command to submit job
       submit: 'qsub'
       # Command to check job status
       stat: 'qstat'
       # Relative path to submit script template
       # The path has to be relative to:
       # $ dirname $(which uap)
       template: 'cluster/submit-scripts/qsub-template.sh'
       # way to define job dependencies
       hold_jid: '-hold_jid'
       # Separator for job dependencies
       hold_jid_separator: ';'
       # Option to set job names
       array_job: '-t 1-%s'
       # Option to submit an array job
       array_job_wquota: '-t 1-%s -tc %s'
       # Options to submit an array job with a quota
       set_job_name: '-N'
       # Option to set path of stderr file
       set_stderr: '-e'
       # Option to set path of stdout file
       set_stdout: '-o'
       # Regex to extract Job ID after submission
       parse_job_id: 'Your job (\d+)'

   # Configuration for a SLURM cluster engine
   slurm:
       identity_test: ['sbatch', '--version']
       identity_answer: 'slurm'
       submit: 'sbatch'
       stat: 'squeue'
       template: 'cluster/submit-scripts/sbatch-template.sh'
       hold_jid: '--dependency=afterany:%s'
       hold_jid_separator: ':'
       array_job: '--array=1-%s'
       array_job_wquota: '--array=1-%s%%%s'
       set_job_name: '--job-name=%s'
       set_stderr: '-e'
       set_stdout: '-o'
       parse_job_id: 'Submitted batch job (\d+)'


Let's browse over the options which need to be set per cluster engine:

``identity_test:``
    Command used to determine if **uap** has been started on a system running
    a cluster engine e.g. ``sbatch --version``.

``identity_answer:``
    **uap** checks if the output of the ``identity_test`` command starts with
    this value e.g. ``slurm``.
    If that is true the cluster type has been detected.

``submit:``
    Command to submit a job onto the cluster e.g. ``sbatch``.

``stat:``
    Command to check the status of jobs on the cluster e.g. ``squeue``.

``template:``
    Path to the submit script template which has to be used for this cluster
    type e.g. ``cluster/submit-scripts/sbatch-template.sh``.


``hold_jid:``
    Option given to the ``submit`` command to define dependencies between
    jobs e.g. ``--dependency=afterany:%s``.
    Placeholder ``%s`` gets replaced with the jobs this job depends on if
    present.

``hold_jid_separator:``
    Separator used to concatenate multiple jobs for ``hold_jid`` e.g. ``:``.

``array_job``:
    Option given to the ``submit`` command to use array jobs e.g.
    ``--array=1-%s``.
    ``%s`` is replaced by the number of jobs.

``array_job_wquota``:
    Option given to the ``submit`` command to use array jobs with quota
    e.g. ``--array=1-%s%%%s`` (will be ``--array=1-100%5`` for *100*
    jobs with a quota of *5*).
    The first ``%s`` is replaced by the number of jobs and the second
    ``%s`` by the quota (if above 0). A literal "%" has to be written
    as ``%%``.

``array_task_id``
    The name of the environment variable set by the resource manager
    that contains the job array id e.g.
    ``SLURM_ARRAY_TASK_ID`` or ``SGE_TASK_ID``.

``set_job_name:``
    Option given to the ``submit`` command to set the job name e.g.
    ``--job-name=%s``.
    ``%s`` is replaced by the job name if present.

``set_stderr:``
    Option given to the ``submit`` command to set the name of the stderr file
    e.g. ``-e``.

``set_stdout:``
    Option given to the ``submit`` command to set the name of the stdout file
    e.g. ``-o``.

``parse_job_id:``
    Python regular expression whose first parenthesized subgroup represents
    the cluster job ID e.g. ``Submitted batch job (\d+)``.

.. _submit_template:

Submit Script Template
======================

The submit script template contains a lot of placeholders which are replaced
if a job is submitted to the cluster with the actual commands.

The submit script templates reside at::

    $ ls $(dirname $(which uap))/cluster/submit-scripts/*
    qsub-template.sh
    sbatch-template.sh

Feel free to add your own templates.
The templates need to contain the following placeholders:

.. _submit_template_submit_options:

``#{SUBMIT_OPTIONS}``
    Will be replaced with the steps ``_cluster_submit_options`` value (see
    :ref:`_cluster_submit_options <_config_file_cluster_submit_options>`), if
    present, or the ``default_submit_options`` value.

.. _submit_template_pre_job_command:

``#{PRE_JOB_COMMAND}``
   Will be replaced with the steps ``_cluster_pre_job_command`` value (see
   :ref:`_cluster_pre_job_command <_config_file_cluster_pre_job_command>`), if
   present, or the ``default_pre_job_command`` value.

.. _submit_template_array_jobs:

``#{ARRAY_JOBS}``
   Will be replaced with a space seperated list of tasks. The resulting array
   will be used in the command for the ``<run ID>`` if the submitted job is
   an array job.

.. _submit_template_command:

``#{COMMAND}``
   Will be replaced with ``uap <project-config>.yaml run-locally <run ID>``.

.. _submit_template_post_job_command:

``#{POST_JOB_COMMAND}``
   Will be replaced with the steps ``_cluster_post_job_command`` value (see
   :ref:`_cluster_post_job_command <_config_file_cluster_post_job_command>`), if
   present, or the ``default_post_job_command`` value.

The submit script template is required by
:ref:`submit-to-cluster <uap-submit-to-cluster>` for job submission to the
cluster.


.. .. [1] |pyyaml_link|

.. |uge_link| raw:: html

   <a href="http://www.univa.com/products/" target="_blank">UGE</a>

.. |slurm_link| raw:: html

   <a href="http://slurm.schedmd.com/" target="_blank">SLURM</a>

.. |yaml_link| raw:: html

   <a href="http://www.yaml.org/" target="_blank">YAML</a>

.. |pyyaml_link| raw:: html

   <a href="http://pyyaml.org/ticket/128" target="_blank">PyYAML</a>

.. |pythex_link| raw:: html

   <a href="http://pythex.org" target="_blank">pythex.org</a>

.. |lmod_link| raw:: html

   <a href="https://lmod.readthedocs.io/en/latest/" target="_blank">lmod</a>

.. |coreutils| raw:: html

    <a href="https://www.gnu.org/software/coreutils/manual/coreutils.html" target="_blank">GNU Core Utilities</a>
