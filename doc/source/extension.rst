..
  This is the documentation for uap. Please keep lines under 80 characters if
  you can and start each sentence on a new line as it decreases maintenance
  and makes diffs more readable.

.. title:: Extending **uap**

..
  This document describes how **uap** can be extended with new analysis steps.


Extending **uap** Functionality
===============================


Implement your own steps
------------------------

The provided pipeline can be easily extended by implementing new steps and
sources. Therefore one does need some basic python programming skills. To add a
new processing step, a single Python file must be placed in ``include/step``
which defines a class with a constructor and two functions. The constructor
(``__init__``) checks for the availability of required tools and tells the
pipeline which connections this step expects (``in/``) and which it provides
(``out/``). The first of the functions  (``setup_runs``) is used for planning all
jobs based on a list of input files or runs and possibly additional information
from previous steps and the second function (``execute``) is used to execute a
specific job. The basic scaffold is shown below.

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

Best practices
**************

There are a couple of things which should be kept in mind when implementing new 
steps or modifying existing steps:

* Make sure errors already show up in ``setup_runs`` instead of ``execute``.
  Therefore look out for things that may fail in ``setup_runs``. Stick to *fail
  early, fail often*. That way errors show up before submitting jobs to the
  cluster and wasting precious cluster waiting time is avoided. 
* Use the ``info`` entry in the returned ``output_run_info`` structure to pass
  information gathered in ``setup_runs`` to ``execute``.
* Likewise, make sure that the tools you'll need in ``execute`` are available.
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



Add the new step to your configuration
--------------------------------------

To insert a new step in a pipeline it has to be added into the ``config.yaml``.
