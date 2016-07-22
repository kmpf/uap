..
  This is the documentation for uap. Please keep lines under 80 characters if
  you can and start each sentence on a new line as it decreases maintenance
  and makes diffs more readable.

.. title:: Extension of uap

..
  This document describes how **uap** can be extended with new analysis steps.

.. _extending-uap:

###############################
Extending **uap** Functionality
###############################

*******************
Implement new steps
*******************

**uap** can be easily extended by implementing new
:ref:`source <config_file_source_steps>` or
:ref:`processing <config_file_processing_steps>` steps.
That requires basic python programming skills.
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
   import <module-you-need>

   # Thirdly import local application files
   from abstract_step import AbstractStep # or AbstractSourceStep

   # Get application wide logger
   logger = getLogger("uap_logger")


Essential imports are the ``from logging import getLogger`` and
``from abstract_step import ...``.
The former is necessary to get access to the application wide logger and
the latter to be able to inherit either from ``AbstractStep`` or
``AbstractSourceStep``.

.. _extending_class_init:

Step 2: Class Definition and Constructor
========================================

Now you need to define a class (which inherits either from ``AbstractStep`` or
``AbstractSourceStep``) and its ``__init__`` method.

.. code-block:: python

   # Either inherit from AbstractSourceStep or AbstractStep
   # class NameOfNewSourceStep(AbstractSourceStep):
   class HelloWorld(AbstractStep):
       # Overwrite constructor
       def __init__(self, pipeline):
           # Call super classes constructor
           super(HelloWorld, self).__init__(pipeline)



The new class needs to be derived from either ``AbstractStep``, for processing
steps, or ``AbstractSourceStep``, for source steps.

The ``__init__`` method should contain the declarations of:

  * Tools used in the step: ``self.require_tool('tool_name')``
  * input and output connection(s): ``self.add_connection('in/*')`` or 
    ``self.add_connection('out/*')``
  * options: ``self.add_option()``

Tools :
  Normally, steps use tools to perform there task.
  Each tool that is going to be used by a step needs to be requested via the
  method ``require_tool('tool_name')``.
  When the step is executed  **uap** searches for ``tool_name`` in the tools
  section of the configuration and uses the information given there to verify
  the tools accessibility.

Connections:
  They are defined by the method ``add_connection(...)``.
  Information is transferred from one step to another via these connections.
  An output connection (``out/something``) of a predecessor step is
  automatically connected with an input connection of the same name
  (``in/something``).
  Connection names should describe the data itself **NOT** the data type.
  For instance, use ``in/genome`` over ``in/fasta``.
  The data type of the input data should be checked in the step anyway to
  execute the correct corresponding commands.

Options:
  Steps can have any number of options.
  Options are defined by the method ``add_option()``.
  There are a bunch of parameters which can be set to specify the option.


.. code-block:: python

   ..
   # Either inherit from AbstractSourceStep or AbstractStep
   # class NameOfNewSourceStep(AbstractSourceStep):
   class NameOfNewProcessingStep(AbstractStep):
       # Overwrite constructor
       def __init__(self, pipeline):
           # Call super classes constructor
           super(NameOfNewProcessingStep, self).__init__(pipeline)

           # Define connections
           self.add_connection('in/some_incoming_data')
           self.add_connection('out/some_outgoing_data')

           # Request tools
           self.require_tool('cat')

           # Add options
           self.add_option('some_option', str, optional=False, 
                           description='Mandatory option')

The single function  ``runs`` is used to plan all jobs based on a list of input
files or runs and possibly additional information from previous steps.
The basic scaffold is shown below.

.. code-block:: python

    import sys
    from abstract_step import *
    import pipeline
    import re
    import process_pool
    import yaml
    
    class Macs14(AbstractStep):
        
        # the constructor
        def __init__(self, pipeline):
            super(Macs14, self).__init__(pipeline)

            # define in and out connections the strings have to start with 'in/'
            # or 'out/'
            self.add_connection('in/something')
            self.add_connection('out/tag1')
            self.add_connection('out/tag2')
            ...
    
            self.require_tool('cat4m')
            self.require_tool('pigz')
            ...

        # all checks of options and input values should be done here
        def setup_runs(self, complete_input_run_info, connection_info):
            # a hash containing information about this step
            output_run_info = {}

            # analyze the complete_input_run_info hash provided by the pipeline
            for step_name, step_input_info in complete_input_run_info.items():
                for input_run_id, input_run_info in step_input_info.items():
                   # assemble your output_run_info
                   # output_run_info has to look like this
                   output_run_info:
                       run_id_1:
                           "output_files":
                               tag1:
                                   output_file_1: [input_file_1, input_file_2, ...]
                                   output_file_2: [input_file_1, input_file_2, ...]
                               tag2:
                                   output_file_3: [input_file_1, input_file_2, ...]
                                   output_file_4: [input_file_1, input_file_2, ...]
                           "info":
                               ...
                           more:
                               ...
                           keys:
                               ...
                       run_id_2:
                           ...

            return output_run_info
        
        # called to actually launch the job (run_info is the hash returned from
        # setup_runs)
        def execute(self, run_id, run_info):
    
            with process_pool.ProcessPool(self) as pool:
                with pool.Pipeline(pool) as pipeline:
                    # assemble the steps pipline here
                    pipeline.append(...)
                    ...
                    # finally launch it
                    pool.launch(...)

The code shown above is the framework for a new step. The most essential part is
the hash returned by setup_runs(), here called ``output_run_info``.

:``run_id``:
    It has to be the unique name of a run (obviously, because its a key value).
    ``output_run_info`` can contain multiple ``run_id`` hashes.

:``"output_files"``:
    This is the only hash key that has to have a fix name. This is used to link
    input to output files.

:``tag[12]``:
    Every ``tag`` has to match ``\w+$`` in the string ``'out/tag'``, which was
    given to ``self.add_connection('out/tag')``. This can be any string, but it
    has to match with the last part of the connection string.

:``output_file_\d``:
    Each ``tag`` has to contain at least one such key. It has to be the name of
    the output file produced by the connection ``'out/tag'``. The value of this
    has to be a list of related input files. The list can have any number of
    entries even zero. Multiple ``output_file_\d`` can rely on the same set of
    input files.

Also very important is to understand the concept of *connections*. They provide
input files prior steps created already. The names of the connections can be
arbitrarily chosen, but should **not** describe the file format but more general
terms. For example an ``out/alignment`` can provide gzipped SAM or BAM files. So
you have to check in setup runs for the file type provided by a connection and
react accordingly. Inspect ``complete_input_run_info`` to find out what your
step gets as input.

.. _uap_tools:

uap tools
=========

You will need to run bash commands like cat, pigz or something else in python. 
In this cases use the uap tool ``exec_group`` (see ``run::new_exec_group()``)

For example you want to separate multiple lines with a specific string out of a
file in a new output file and in addition to this copy the output file.
A possible bash way is:

.. code-block:: bash

    $ cat source_file | grep search_string > output_file
    $ cp output_file new_file

For sure, for this task grep would be sufficient. But for the example we want to use a pipe.

Now the uap way:

.. code-block:: python

    # create an new exec_group object
    exec_group = run.new_exec_group()

    # create an output file for the pipeline
    cat_out = run.add_output_file(
        'file',
        '%s.txt' % (run_id),
        [input_path])

    # create a command with cat and grep combined through pipe
    with exec_group.add_pipeline() as cat_pipe:
        # create the cat command
        cat_command = [self.get_tool('cat'), input_path]

        # create the grep command
        search_string = 'foobar'
        grep_command = [self.get_tool('grep'), search_string]

        # add commands to the command pipeline
        cat_pipe.add_command(cat_command)
        cat_pipe.add_command(grep_command, stdout_path= cat_out)

    # create a copy output file
    cp_out = run.add_output_file(
        'file',
        '%s_copy.txt' % (run_id),
        [input_path])

    # create copy command
    cp_command = [self.get_tool('cp'), cat_out, cp_out]

    # add copy command to the pipeline
    exec_group.add_command(cp_command)

All the single commands will be collected and uap will execute the command list in the specified order.

Best practices
==============

There are a couple of things which should be kept in mind when implementing new 
steps or modifying existing steps:

* Make sure errors already show up in ``runs``.
  So, look out for things that may fail in ``runs``.
  Stick to *fail early, fail often*.
  That way errors show up before submitting jobs to the cluster and wasting 
  precious cluster waiting time is avoided. 
* Make sure that the tools you'll need in ``runs`` are available.
  Check for the availability of tools within the constructor ``__init__``.

.. code-block:: python

    # make sure tools are available
    self.require_tool('pigz')
    self.require_tool('cutadapt')

* Make sure your disk access is as cluster-friendly as possible (which 
  primarily means using large block sizes and preferably no seek operations). 
  If possible, use ``unix_pipeline`` to wrap your commands in ``pigz``, ``dd``,
  or ``cat4m`` with a large block size like 4 MB. 
  Although this is not possible in every case (for example when seeking 
  in files is involved), it is straightforward with tools that read a 
  continuous stream from ``stdin`` and write a continuous stream to 
  ``stdout``.
* **NEVER**  remove files! If files need to be removed report the issue and 
  exit **uap**. Only the user should delete files.
* Always use ``os.path.join(...)`` when you handle paths.
* Use bash commands like ``mkfifo`` over python library equivalents like
  ``os.mkfifo()``
* If you need to decide between possible ways to implement a step, stick to the
  more flexibel (often more configuration extensive one).
  You don't know what other user might need, so let them decide.

**************************************
Add the new step to your configuration
**************************************

To make a new step known to **uap**, it has to be copied into either of these
folders:

``include/sources/``
  for all source steps

``include/steps/``
  for all processing steps

If the Python step file exist at the correct location the step needs to be added
to the YAML configuration file as described in :doc:`configuration`.

