..
  This is the documentation for uap. Please keep lines under
  80 characters if you can and start each sentence on a new line as it 
  decreases maintenance and makes diffs more readable.

.. title:: Configuration of uap

..
  This document aims to describe how to configure **uap**.

.. _configuration-of-uap:

***************************
Analysis Configuration File
***************************

**uap** operates on |yaml_link| files which define data
analysis.
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

``destination_path``
--------------------

The value of ``destination_path`` is the directory where **uap** is going
to store the created files.

.. It is possible to use a different directory for volatile files (see ).

.. code-block:: yaml

    destination_path: "/path/to/workflow/output"


``constants``
-------------

This section is the place where repeatedly used constants should be defined.
For instance absolute paths to the genome index files can be defined as
constant.

.. code-block:: yaml

   - &genome_faidx
        genomes/animalia/chordata/mammalia/primates/homo_sapiens/hg19/samtools_faidx/hg19_all_chr_UCSC-download-B7ceRp9K/hg19_all_chr_UCSC-download.fasta.fai

Later on the value can be reused by typing ``*genome_faidx``.
**There are no restrictions about what can be defined here.**

.. _config-file-steps:

``steps``
---------

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
an **``uap``** analaysis.

.. _config_file_source_steps:

Source Steps
^^^^^^^^^^^^
They provide input files for the analysis.
They might start processes such as downloading files or demultiplexing
sequence reads.
But, they do not have dependencies, they can introduce files from outside the
destination path (see `Destination_path Section`_), and they are
usually the first steps of an analysis.

For example if you want to work with fastq files, the first step is to import the required files. For this task the source step fastq_source is the right solution.

A possible step definition could look like this:

.. code-block:: yaml

    steps:
        input_step (fastq_source):
        pattern: /Path/to/fastq-files/*.gz
        group: ([SL]\w+)_R[12]-00[12].fastq.gz
        sample_id_prefix: MyPrefix
        first_read: '_R1'
        second_read: '_R2'
        paired_end: True

The single keys will be described at :doc:`steps`. For defining the ``group`` key a regular expression is used. If you are not familiar with this you can read about it and test your regular expression at |pythex_link|.

.. _config_file_processing_steps:

Processing Steps
^^^^^^^^^^^^^^^^

They depend upon one or more predecessor steps and work with their output
files.
Output files of processing steps are automatically named and placed by **uap**.
Processing steps are usually configurable.
For a complete list of available options please visit :doc:`steps` or use the
subcommand :ref:`uap-steps`.

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
  of its succeeding step data gets passed on automatically.
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
            _depends:
                - chr1
                - chr2
            # Equivalent to:
            # _depends: [chr1, chr2]
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
The message is output if the :ref:`status <uap-status>` is checked and looks like this::

  Hint: You could save 156.9 GB of disk space by volatilizing 104 output files.
  Call 'uap <project-config>.yaml volatilize --srsly' to purge the files.

If the user executes the :ref:`volatilize <uap-volatilize>` command the output
files are replaced by placeholder files.

.. _config_file_cluster_submit_options:

**_cluster_submit_options**

    This string contains the entire submit options which will be set in the
    submit script.
    This option allows to overwrite the values set in 
    :ref:`config_file_default_submit_options`.

.. _config_file_cluster_pre_job_command:

**_cluster_pre_job_command**

    This string contains command(s) that are executed **BEFORE uap** is started
    on the cluster.
    This option allows to overwrite the values set in 
    :ref:`config_file_default_pre_job_command`.

.. _config_file_cluster_post_job_command:

**_cluster_post_job_command**

    This string contains command(s) that are executed **AFTER uap** did finish
    on the cluster.
    This option allows to overwrite the values set in 
    :ref:`config_file_default_post_job_command`.

.. _config_file_cluster_job_quota:

**_cluster_job_quota**

    This positive number defines the number of jobs of the same type that can
    run simultaneously on a cluster.
    This option allows to overwrite the values set in 
    :ref:`config_file_default_job_quota`.

.. _uap_config_tools:

``tools``
---------

The ``tools`` section must list all programs required for the execution of a
particular analysis.
**uap** uses the information given here to check if a tool is available given
the current environment.
This is particularly useful on cluster systems were software might not always
be loaded.
Also, **uap** logs the version of each tool used by a step.

By default, version determination is simply attempted by calling the program
without command-line arguments.

If a certain argument is required, specify it in ``get_version``. 

If a tool does not exit with code 0, you can find out which code is it.
Execute the required command and after this type ``echo $?`` in the same shell.
The output is the exit code of the last executed command.
You can use it to specify the exit code in ``exit_code``.

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

If you are working on a cluster running |uge_link|
or |slurm_link| you can also use their module system.
You need to know what actually happens when you load or unload a module::

  $ module load <module-name>
  $ module unload <module-name>

As far as I know is ``module`` neither a command nor an alias.
It is a BASH function. So use ``declare -f`` to find out what it is actually
doing::

  $ declare -f module

The output should look like this:

.. code-block:: bash

    module ()
        {
            eval `/usr/local/modules/3.2.10-1/Modules/$MODULE_VERSION/bin/modulecmd bash $*`
        }

An other possible output is:

.. code-block:: bash

    module () 
        { 
            eval $($LMOD_CMD bash "$@");
            [ $? = 0 ] && eval $(${LMOD_SETTARG_CMD:-:} -s sh)
        }

In this case you have to look in ``$LMOD_CMD`` for the required path::

    $ echo $LMOD_CMD

Now you can use this newly gathered information to load a module before use
and unload it afterwards.
You only need to replace ``$MODULE_VERSION`` with the current version of the
module system you are using and ``bash`` with ``python``.
A potential ``bedtools`` entry in the ``tools`` section, might look like this.

.. code-block:: yaml

    tools:
        ....
        bedtools:
            module_load: /usr/local/modules/3.2.10-1/Modules/3.2.10/bin/modulecmd python load bedtools/2.24.0-1
            module_unload: /usr/local/modules/3.2.10-1/Modules/3.2.10/bin/modulecmd python unload bedtools/2.24.0-1
            path: bedtools
            get_version: --version
            exit_code: 0


.. NOTE:: Use ``python`` instead of ``bash`` for loading modules via **uap**.
          Because the module is loaded from within a python environment and
          not within a BASH shell.

.. _config_file_cluster: 

``cluster``
-----------

The value of ``cluster`` is needed if the analysis is executed on a cluster,

.. code-block:: yaml

    cluster:
        default_submit_options: "-pe smp #{CORES} -cwd -S /bin/bash -m as -M me@example.com -l h_rt=1:00:00 -l h_vmem=2G"
        default_pre_job_command: "echo 'Started the run!'"
        default_post_job_command: "echo 'Finished the run!'"
        default_job_quota: 5

.. _config_file_default_submit_options:

**default_submit_options**

.. _config_file_default_pre_job_command:

**default_pre_job_command**

.. _config_file_default_post_job_command:

**default_post_job_command**

.. _config_file_default_job_quota:

**default_job_quota:**

Example Configurations
======================

Please check out the example configurations provided inside the ``example-configurations`` folder of **uap**'s installation directory.


**************************
Cluster Configuration File
**************************

The cluster configuration file resides at::

    $ ls -la $(dirname $(which uap))/cluster/cluster-specific-commands.yaml

This YAML file contains a dictionary per cluster type, that looks like that::

    uge: # Uniq name of the cluster engine
        identity_test: ['qstat', '-help'] # Command to get version information
        identity_answer: 'UGE' # The output of the above command for that cluster
        submit: 'qsub' # Command to submit job
        stat: 'qstat' # Command to check job status
        template: 'cluster/submit-scripts/qsub-template.sh' # Path to template for submit script (relative to dirname $(which uap))
        hold_jid: '-hold_jid' # way to define job dependencies
        hold_jid_separator: ';' # Separator for job dependencies
        set_job_name: '-N' # Way to set job names
        set_stderr: '-e' # Way to set path to file for stderr
        set_stdout: '-o' # Way to set path to file for stdout
        parse_job_id: 'Your job (\d+)' # Regex to extract Job ID after submission


Ausbauen!!!




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
