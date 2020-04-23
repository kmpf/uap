..
  This is the documentation for uap. Please keep lines under
  80 characters if you can and start each sentence on a new line as it 
  decreases maintenance and makes diffs more readable.

.. title:: Information for uap Developers

..
  This document describe different issues concerning the development of **uap**.


******************************
Information for uap Developers
******************************

Documentation
=============

The official documentation is hosted on readthedocs (add link). But as
developer you need to be able to create the documentation locally. So let's
focus at first on:

Create Documentation Locally
----------------------------

For testing purposes its necessary to be able to build the documentation on your
local machine. Follow the instructions below.

Prerequisites
~~~~~~~~~~~~~

Before the documentation can be build, we need to install some packages and
libraries.

The documentation is build with Sphinx, so install it (for Ubuntu)::

  $ sudo aptitude install python-sphinx

The uap documentation requires the ``psutil`` python library to be available
in your global python environment::

  $ sudo pip install psutil

It also requires the readthedocs theme installed (for Ubuntu)::

  $ sudo aptitude install python-sphinx-rtd-theme sphinx-rtd-theme-common


Build Documentation
~~~~~~~~~~~~~~~~~~~

The documentation can be build by::

  $ cd doc
  $ make html

This should do the trick.

Create Documentation on RTD
---------------------------

Commit your changes to one of the branches that are configured for automatic
builds via RTD. Push your commits to the yigbt-repository. This should do the
trick.

If you are not allowed to directly push your commits. Please provide a pull
request via github. This should do the trick as well.

Continuous Integration
======================

At the moment we use Travis CI as continuous integration platform. The goal is
to test as many step use cases as possible. Travis installs uap in a clean
Docker image. Installation instructions for software that needs to be installed
additionally, needs to be provided  you want to install 

.travis.yml
-----------

This file contains the test matrix, the commands executed before installing uap,
the commands executed to install uap and the commands which are executed to
perform the tests you are interested in.

At the moment most of the subcommands are tested on Travis CI using a dummy
workflow.

travis_uap_config.yaml
~~~~~~~~~~~~~~~~~~~~~~

The dummy workflow we use to test steps is defined in the file
``example-configurations/travis-example/travis_uap_config.yaml``.
So if you create a new step, please add it to this configuration file. Provide
as many step configurations as you think are appropriate.

