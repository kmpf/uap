..
  This is the documentation for uap. Please keep lines under
  80 characters if you can and start each sentence on a new line as it 
  decreases maintenance and makes diffs more readable.

.. title:: Configuration of uap

..
  This document aims to describe how to configure **uap**.

.. _configuration-of-uap:
##################
Configuration File
##################

**uap** operates on |yaml_link| files which define data
analysis.
These files are called configuration files.

A configuration file describes a analysis completely.
Configurations consist of four sections (let's just call them sections,
although technically, they are keys):

* ``destination_path`` -- points to the directory where the result files,
  annotations and temporary files are written to
* ``email`` -- when submitting jobs on a cluster, messages will be sent to 
  this email address by the cluster engine (nobody@example.com by default)
* ``constants`` -- defines constants for later use (define repeatedly used
  values as constants to increase readability of the following sections)
* ``steps`` -- defines the source and processing steps and their order 
* ``tools`` -- defines all tools used in the analysis and how to determine 
  their versions (for later reference)

If you want to know more about the notation that is used in this file, have a
closer look at the |yaml_link| definition.

********************************
Sections of a Configuration File
********************************

.. _config-file-destination-path:
Destination_path Section
========================

The value of ``destination_path`` is the directory where **uap** is going
to store the created files.

.. It is possible to use a different directory for volatile files (see ).

.. code-block:: yaml

    destination_path: "/path/to/uap/output"

Email Section
=============

The value of ``email`` is needed if the analysis is executed on a cluster,
which can use it to inform the person who started **uap** about status
changes of submitted jobs.

.. code-block:: yaml

    email: "your.name@mail.de"


Steps Section
=============

The ``steps`` section is the core of the analysis file, because it defines when
steps are executed and how they depend on each other.
All available steps are described in detail in the steps documentation: 
:doc:`steps`.
This section (technically it is a dictionary) contains a key for every step,
therefore each step must have a unique name [1]_.
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
            
There are two different types of steps:

.. _config_file_source_steps:
Source Steps
------------
They provide input files for the analysis.
They might start processes such as downloading files or demultiplexing
sequence reads.
But, they do not have dependencies, they can introduce files from outside the
destination path (see :ref:`config_file_destination_path`), and they are
usually the first steps of an analysis.

.. _config_file_processing_steps:
Processing Steps
----------------
They depend upon one or more predecessor steps and work with their output
files.
Output files of processing steps are automatically named and placed by **uap**.
Processing steps are usually configurable.
For a complete list of available options please visit :doc:`steps` or use the
subcommand :ref:`uap_steps`.

Reserved Keywords for Steps
---------------------------

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


.. _uap-volatile:
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

.. _uap_config_tools_section:
Tools Section
=============

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
If a tool does not exit with exit code 0, find out which code it is by typing
``echo $?`` into Bash and specify the exit code in ``exit_code``.

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

**********************
Example Configurations
**********************

Please check out the example configurations provided inside the ``example-configurations`` folder of **uap**'s installation directory.

.. [1] |pyyaml_link|

.. |uge_link| raw:: html

   <a href="http://www.univa.com/products/" target="_blank">UGE</a>.

.. |slurm_link| raw:: html

   <a href="http://slurm.schedmd.com/" target="_blank">SLURM</a>.

.. |yaml_link| raw:: html

   <a href="http://www.yaml.org/" target="_blank">YAML</a>.

.. |pyyaml_link| raw:: html

   <a href="http://pyyaml.org/ticket/128" target="_blank">PyYAML does not complain about duplicate keys</a>.
