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

The installation requires |virtual_env_link|, |git_link| and
|graphviz_link|.
So, please install it if its not already installed.::

  $ sudo apt-get install python-virtualenv git graphviz

**uap** does **NOT** include any tools necessary for the data analysis.
It is expected that the required tools are **already installed**.

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
Afterwards it installs |py_yaml_link|, |num_py_link|, |bio_python_link| and
|psutil_link| into the freshly created environment.
There is no harm in accidentally running this script multiple times.

*********************************
Making **uap** Globally Available
*********************************

**uap** can be used globally.
On Unix-type operating systems it is advised to add the installation path to
your ``$PATH`` variable.
Therefore change into the **uap** directory and execute::

.. code-block:: bash

  $ echo ""PATH=$PATH:$(pwd)" >> ~/.bashrc 
  $ source ~/.bashrc
  OR
  $ echo ""PATH=$PATH:$(pwd)" >> ~/.bash_profile
  $ source ~/.bash_profile


.. |github_uap_link| raw:: html

   <a href="https://github.com/kmpf/uap" target="_blank">uap's github repository</a>.

.. |virtual_env_link| raw:: html

   <a href="https://pypi.python.org/pypi/virtualenv" target="_blank">virtualenv</a>

.. |git_link| raw:: html

   <a href="https://git-scm.com/" target="_blank">git</a>

.. |github_uap_link| raw:: html

   <a href="https://github.com/kmpf/uap" target="_blank">uap's github repository</a>.

.. |graphviz_link| raw:: html

   <a href="http://www.graphviz.org/" target="_blank">graphviz</a>

.. |py_yaml_link| raw:: html
 
    <a href="https://pypi.python.org/pypi/PyYAML" target="_blank">PyYAML</a>

.. |num_py_link| raw:: html
 
    <a href="https://pypi.python.org/pypi/numpy" target="_blank">NumPy</a>

.. |bio_python_link| raw:: html
 
    <a href="https://pypi.python.org/pypi/biopython" target="_blank">biopython</a>

.. |psutil_link| raw:: html
 
    <a href="https://pypi.python.org/pypi/psutil" target="_blank">psutil</a>
