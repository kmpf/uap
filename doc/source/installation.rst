..
  This is the documentation for rnaseq-pipeline. Please keep lines under
  80 characters if you can and start each sentence on a new line as it 
  decreases maintenance and makes diffs more readable.

.. title:: Installation of uap

..
  This document aims to describe how to install **uap**.

.. _installation-of-uap:
#######################
Installation of **uap**
#######################

*************
Prerequisites
*************

The installation requires Python's ``virtualenv`` tool.
So, please install it if its not already installed.::

  $ sudo apt-get install python-virtualenv
  OR
  $ sudo pip install virtualenv

************************
Downloading the Software
************************

Download the software from |github_uap_link| like this::

  $ git clone https://github.com/kmpf/uap.git

*****************************    
Setting Up Python Environment
*****************************

After cloning the repository, change into the created directory and run the 
bootstrapping script ``bootstrap.sh``::

  $ cd uap
  $ ./bootstrap.sh

The script creates the required Python environment (which will be located in
``./python_env/``).
Afterwards it installs
`PyYAML <https://pypi.python.org/pypi/PyYAML>`_,
`NumPy <https://pypi.python.org/pypi/numpy>`_,
`biopython <https://pypi.python.org/pypi/biopython>`_, and
`psutil <https://pypi.python.org/pypi/psutil>`_ into the freshly created
environment.
There is no harm in accidentally running this script multiple times.

*********************************
Making **uap** Globally Available
*********************************

**uap** can be used globally.
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

.. |github_uap_link| raw:: html

   <a href="https://github.com/kmpf/uap" target="_blank">uap's github repository</a>.

.. |github_uap_link| raw:: html

   <a href="https://github.com/kmpf/uap" target="_blank">uap's github repository</a>.
.. |github_uap_link| raw:: html

   <a href="https://github.com/kmpf/uap" target="_blank">uap's github repository</a>.
.. |github_uap_link| raw:: html

   <a href="https://github.com/kmpf/uap" target="_blank">uap's github repository</a>.
.. |github_uap_link| raw:: html

   <a href="https://github.com/kmpf/uap" target="_blank">uap's github repository</a>.
