..
  This is the documentation for rnaseq-pipeline. Please keep lines under
  80 characters if you can and start each sentence on a new line as it 
  decreases maintenance and makes diffs more readable.

.. title:: Installing **uap**

..
  This document aims to describe how to install **uap**.

Installation
============


Downloading the software
------------------------

The repository can be obtained like this::

    $ git clone git@github.com:kmpf/uap.git
    
We'll need ``virtualenv``, so install it if you don't have it already::

    $ sudo apt-get install python-virtualenv

After cloning the repository, change into the created directory and run the 
bootstrapping script to create the required Python environment (which will be
located in ``./python_env/``)::

    $ cd uap
    $ ./bootstrap.sh

There's no harm in accidentally running this script multiple times. 
Also, it will compile ``cat4m``, a tool which can be found at 
``./tools/cat4m`` and which is able to read arbitrary input files in chunks 
of 4 MB and print them to stdout (we'll need this often in the pipeline,
as ``cat`` reads in system-default blocks of 32 kB which is ok for a normal
system but leads to high I/O load on a cluster system).

Making **uap** Globally Available
---------------------------------

You can use **uap** system-wide.
On Unix-type operating systems it is advised to add the installation path to
your ``$PATH`` variable.
Therefore change into the **uap** directory and execute::

    $ pwd

Copy the output and add a line to your ``~/.bashrc`` or ``~/.bash_profile``
like the following and replace ``<uap-path>`` with the copied output:

.. code-block:: bash

    PATH=$PATH:<uap-path>


Finally, make the changes known to your environment by sourcing the changed
file::
    $ source ~/.bashrc
    $ source ~/.bash_profile
