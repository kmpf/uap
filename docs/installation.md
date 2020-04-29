<!--
  This is the documentation for rnaseq-pipeline. Please keep lines under
  80 characters if you can and start each sentence on a new line as it 
  decreases maintenance and makes diffs more readable.
-->

<!--
  This document aims to describe how to install **uap**.
-->

# Installation

## Prerequisites

The installation requires [virtualenv](https://pypi.python.org/pypi/virtualenv),
[git](https://git-scm.com/) and [graphviz](http://www.graphviz.org/).
So, please install it if its not already installed.

```bash
  
$ sudo apt-get install python-virtualenv git graphviz

```

**uap** does **NOT** include any tools necessary for the data analysis.
It is expected that the required tools are **already installed**.

## Clone the Repository

Clone the repository from https://github.com/kmpf/uap like this:

```bash

$ git clone https://github.com/kmpf/uap.git

```

## Set Up Python Environment

Change into the created directory and execute the script `bootstrap.sh`:

```bash

$ cd uap
$ ./bootstrap.sh

```

The script creates the required Python environment (which will be located in
`./python_env/`).
Afterwards it installs [PyYAML](https://pypi.python.org/pypi/PyYAML),
[NumPy](https://pypi.python.org/pypi/numpy),
[biopython](https://pypi.python.org/pypi/biopython)
and [psutil](https://pypi.python.org/pypi/psutil)|psutil_link| into the freshly
created environment.
There is no harm in accidentally running this script multiple times.

## Make **uap** Globally Available

**uap** can be used globally.
On Unix-type operating systems it is advised to add the installation path to
your ``$PATH`` variable.
Therefore change into the **uap** directory and execute::

  $ echo ""PATH=$PATH:$(pwd)" >> ~/.bashrc 
  $ source ~/.bashrc
  OR
  $ echo ""PATH=$PATH:$(pwd)" >> ~/.bash_profile
  $ source ~/.bash_profile




