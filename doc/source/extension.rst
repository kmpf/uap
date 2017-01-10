..
  This is the documentation for uap. Please keep lines under 80 characters if
  you can and start each sentence on a new line as it decreases maintenance
  and makes diffs more readable.

.. title:: Extension of uap

..
  This document describes how **uap** can be extended with new analysis steps.

.. _extending-uap:

#####################
Add New Functionality
#####################

*******************
Implement New Steps
*******************

**uap** can be easily extended by implementing new
:ref:`source <config_file_source_steps>` or
:ref:`processing <config_file_processing_steps>` steps.
This requires basic python programming skills.
New steps are added to **uap** by placing a single Python file into one of these
folders in the **uap** installation directory:

``include/sources``
  Place source step files here

``include/steps``
  Place processing step files here

Let's talk about how to implement such **uap** steps.

.. _extending_import:

Step 1: Import Statements and Logger
====================================

At the beginning of every step please import the required modules and create a
logger object.

.. code-block:: python

   # First import standard libraries
   import os
   from logging import getLogger

   # Secondly import third party libraries
   import yaml

   # Thirdly import local application files
   from abstract_step import AbstractStep # or AbstractSourceStep

   # Get application wide logger
   logger = getLogger("uap_logger")


Essential imports are the ``from logging import getLogger`` and
``from abstract_step import ...``.
The former is necessary to get access to the application wide logger and
the latter to be able to inherit either from ``AbstractStep`` or
``AbstractSourceStep``.

.. _extending_class_def:

Step 2: Class Definition
========================

Now you need to define a class (which inherits either from ``AbstractStep`` or
``AbstractSourceStep``) and its ``__init__`` method.

.. code-block:: python

   class ConcatenateFiles(AbstractStep):
       # Overwrite initialisation
       def __init__(self, pipeline):
           # Call super classes initialisation
           super(ConcatenateFiles, self).__init__(pipeline)

   ..

The new class needs to be derived from either ``AbstractStep``, for processing
steps, or ``AbstractSourceStep``, for source steps.

.. _extending_class_init:

Step 3: ``__init__`` Method
===========================

The ``__init__`` method is the place where you should declare:

Tools via ``self.require_tool('tool_name')``:
  Steps usually require tools to perform their task.
  Each tool that is going to be used by a step needs to be requested via the
  method ``require_tool('tool_name')``.
  **uap** tests the existence of the required tools whenever it constructs the
  directed acyclic graph (DAG) of the analysis.
  The test is based on the information provided in the
  :ref:`tools section <uap_config_tools>` of the
  :ref:`analysis configuration <analysis_configuration>`.
  An entry for ``tool_name`` has to exist and to provide information to verify
  the tools accessibility.

Connections via ``add_connection(...)``:
  Connections are defined by the method ``add_connection(...)``.
  They are used to transfer data from one step to another.
  If a step defines an output connection ``out/something`` and a subsequent
  step defines an input connection named ``in/something``, then the files
  beloging to ``out/something`` will be available via the connection
  ``in/something``.

  Please name connection in a way that they describe the data itself and
  **NOT** the data type.
  For instance, use ``in/genome`` over ``in/fasta``.
  The data type of the received input data should be checked by the steps
  to make sure to execute the correct commands.

  **TODO**: Reanimate the constraints feature. It would often save some lines
  of code to be able to define constraints on the connections.

Options via ``self.add_option()``:
  Options allow to influence the commands executed by a step.
  It is advisable to provide as many meaningful options as possible to keep
  steps flexible.
  Steps can have any number of options.
  Options are defined via the method ``add_option()``.
  
  The ``add_option()`` method allows to specify various information about
  the option.
  The method parameters are these:

  1. ``key``
         name of the option (if possible include the name of the tool
         this option influences e.g. ``dd-blocksize`` to set ``dd`` blocksize)

  2. ``option_type``
         The option type has to be at least one of ``int``, ``float``, ``str``,
         ``bool``, ``list``, or ``dict``.

  3. ``optional`` (Boolean)
         Defines if the option is mandatory (``False``) or optional (``True``).

  4. ``choices``
         List of valid values for the option.

  5. ``default``
         Defines the default value for the option.

  6. ``description``
         The description of the functionality of the option.
         


.. code-block:: python

   ..

           # Define connections
           self.add_connection('in/text')
           self.add_connection('out/text')

           # Request tools
           self.require_tool('cat')

           # Options for workflow
           self.add_option('concatenate_all_files', bool, optional=False,
                           default=False, description="Concatenate all files from "
                           "all runs, if 'True'.")

           # Options for 'cat' (see manpage)
           self.add_option('show-all', bool, optional=True,
                           description="Show all characters")
                           
           self.add_option('number-nonblank', int, optional=True,
                           description="number nonempty output lines, "
                           "overrides --number")

           self.add_option('show-ends', bool, optional=True,
                           description="display $ at end of each line")

           self.add_option("number", int, optional=True,
                           description="number all output lines")

           self.add_option("squeeze-blank", bool, optional=True,
                           description="suppress repeated empty output lines")

           self.add_option("show-tabs", bool, optional=True,
                           description="display TAB characters as ^I")

           self.add_option("show-nonprinting", bool, optional=True,
                            description="use ^ and M- notation, except for "
                            "LFD and TAB")

   ..

.. _extending_class_runs:

Step 4: ``runs`` Method
=======================

The ``runs`` method is where all the work is done.
This method gets handed over a dictionary of dictionaries.
The keys of the first dictionary are the run IDs (often resembling the samples).
The values of the first dictionary is another dictionary.
The keys of that second dictionary are the connections e.g. "in/text" and the
values are the corresponding files belonging to that connection.

Let's inspect all the run IDs, connections, and input files we got from our
upstream steps.
And let's tore all files we received in a list for later use.

.. code-block:: python

   ..

       def runs(self, run_ids_connections_files):
           all_files = list()
           # Let's inspect the run_ids_connections_files data structure
           for run_id in run_ids_connections_files.keys():
               logger.info("Run ID: %s" % run_id)
               for connection in run_ids_connections_files[run_id].keys():
                   logger.info("Connection: %s" % connection)
                   for in_file in run_ids_connections_files[run_id][connection]:
                       logger.info("Input file: %s" % in_file)
                       # Collect all files
                       all_files.append(in_file)
   
   ..

It comes in handy to assemble a list with all options for ``cat`` here.

.. code-block:: python

   ..

        # List with options for 'cat'
        cat_options = ['show-all', 'number-nonblank', 'show-ends', 'number',
                       'squeeze-blank', 'show-tabs', 'show-nonprinting']

        # Get all options which were set
        set_options = [option for option in cat_options if \
                       self.is_option_set_in_config(option)]

        # Compile the list of options
        cat_option_list = list()
        for option in set_options:
            # bool options look different than ...
            if isinstance(self.get_option(option), bool):
                if self.get_option(option):
                    cat_option_list.append('--%s' % option)
            # ... the rest ...
            else:
                cat_option_list.append('--%s' % option)
                # ... make sure to cast the values to string
                cat_option_list.append(str(self.get_option(option)))
                
   ..

What should happen if we are told to concatenate all files from all input runs?
We have to create a single run with a new run ID 'all_files'.
The run consists of a ``exec_group`` that runs the ``cat`` command.

.. note::

   An ``exec_group`` is a list of commands which are executed in one go.
   You might create multiple ``exec_group``'s if you need to make sure a set of
   commands finished before another set is started.
   An ``exec_group`` can contain commands and pipelines.
   They can be added like this:

   .. code-block:: python
                   
      # Add a single command
      exec_group.add_command(...)

      # Add a pipeline to an exec_group
      with exec_group.add_pipeline as pipe:
         ...
         # Add a command to a pipeline
         pipe.add_command(...)

The result of the concatenation is written to an output file.
The run object needs to know about each output file that is going to be created.

.. note::

   An output file is announced via the run objects
   ``add_output_file(tag, out_path, in_paths)`` method.
   The method parameters are:

   1. ``tag``: The name of the out connection e.g. 'text' for 'out/text'
   2. ``out_path``: The name of the output file (best practice is to add the
      run ID to the file name)
   3. ``in_paths``: The input files this output file is based on

.. code-block:: python

   ..

        # Okay let's concatenate all files we get
        if self.get_option('concatenate_all_files'):
            run_id = 'all_files'

            # New run named 'all_files' is created here
            with self.declare_run(run_id) as run:

                # Create an exec
                with run.new_exec_group() as exec_group:
                    # Assemble the cat command
                    cat = [ self.get_tool('cat') ]
                    # Add the options to the command
                    cat.extend( cat_option_list )
                    cat.extend( all_files )
                    
                    # Now add the command to the execution group
                    exec_group.add_command(
                        cat,
                        stdout_path = run.add_output_file(
                            'text',
                            "%s_concatenated.txt" % run_id,
                            all_files)
                    )

   ..

What should happen if all files of an input run have to be concatenated?
We create a new run for each input run and concatenate all files that
belong to the input run.

.. code-block:: python

        # Concatenate all files from a runs 'in/text' connection
        else:
            # iterate over all run IDs ...
            for run_id in run_ids_connections_files.keys():
                input_paths = run_ids_connections_files[run_id]['in/text']
                # ... and declare a new run for each of them.
                with self.declare_run(run_id) as run:
                    with run.new_exec_group() as exec_group:
                        # Assemble the cat command
                        cat = [ self.get_tool('cat') ]
                        # Add the options to the command
                        cat.extend( cat_option_list )
                        cat.extend( input_paths )
                        
                        # Now add the command to the execution group
                        exec_group.add_command(
                            cat,
                            stdout_path = run.add_output_file(
                                'text',
                                "%s_concatenated.txt" % run_id,
                                input_paths)
                        )

That's it.
You created your first **uap** processing step.


Step 5: Add the new step to **uap**
===================================

You have to make the new step known to **uap**.
Save the complete file into **uap**'s ``include/steps`` folder.
Processing step files are located at **uap**'s ``include/steps/`` folder
and source step files at **uap**'s ``include/sources/`` folder.

You can control that your step is correctly "installed" if its included in the
list of all source and processing steps::

  $ ls -la $(dirname $(which uap))/include/sources
  ... Lists all available source step files

  $ ls -la $(dirname $(which uap))/include/steps
  ... Lists all available processing step files

You can also use **uap**'s :ref:`steps <uap-steps>` subcommand to get
information about installed steps.

If the step file exists at the correct location that step can be used
in an :ref:`analysis configuration file <analysis_configuration>`.

A potential example YAML file named ``test.yaml`` could look like this:

.. code-block:: yaml

    destination_path: example-out/test/
    
    steps:
        ##################
        ## Source steps ##
        ##################
    
        raw_file_source:
            pattern: example-data/text-files/*.txt
            group: (.*).txt
    
        ######################
        ## Processing steps ##
        ######################
    
        cat:
            _depends: raw_file_source
            _connect:
                in/text:
                    - raw_file_source/raw
            concatenate_all_files: False
    
    tools:
        cat:
            path: cat
            get_version: '--version'
            exit_code: 0

You need to create the destination path and some text files matching the
pattern ``example-data/text-files/*.txt``.
Also you see the work of the ``_connect`` keyword in play.
Check the status of the configured analysis::

  $ uap test.yaml status
  Ready runs
  ----------
  [r] cat/Hello_america
  [r] cat/Hello_asia
  [r] cat/Hello_europe
  [r] cat/Hello_world
  
  runs: 4 total, 4 ready



.. _extending_best_practices:

**************
Best Practices
**************

There are a couple of things you should keep in mind while implementing new 
steps or modifying existing ones:

* **NEVER**  remove files!
  If files need to be removed report the issue and exit **uap** or force the
  user to call a specific subcommand.
  Never delete files without permission by the user.
* Make sure errors already show up in when the steps ``runs()`` method is
  called the first time.
  So, look out for things that may fail in ``runs``.
  Stick to *fail early, fail often*.
  That way errors show up before submitting jobs to the cluster and wasting 
  precious cluster waiting time is avoided.
* Make sure that all tools which you request inside the ``runs()`` method
  are also required by the step via ``self.require_tool()``.
  Use the ``__init__()`` method to request tools.
* Make sure your disk access is as cluster-friendly as possible (which 
  primarily means using large block sizes and preferably no seek operations). 
  If possible, use pipelines to wrap your commands in ``pigz`` or ``dd``
  commands.
  Make the used block size configurable. 
  Although this is not possible in every case (for example when seeking 
  in files is involved), it is straightforward with tools that read a 
  continuous stream from ``stdin`` and write a continuous stream to 
  ``stdout``.
* Always use ``os.path.join(...)`` to handle paths.
* Use bash commands like ``mkfifo`` over python library equivalents like
  ``os.mkfifo()``.
  The ``mkfifo`` command is hashed while an ``os.mkfifo()`` is not.
* Keep your steps as flexible as possible.
  You don't know what other user might need, so let them decide.


Usage of ``dd`` and ``mkfifo``
==============================

**uap** relies often on ``dd`` and FIFOs to process data with fewer
disk read-write operations.
Please provide a step option to adjust the ``dd`` blocksize (this option
is usually called ``dd-blocksize``).
Create your steps in a way that they perform the least filesystem operations.
Some systems might be very sensitive to huge numbers of read-write operations.
