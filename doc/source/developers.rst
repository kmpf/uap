..
  This is the documentation for uap. Please keep lines under
  80 characters if you can and start each sentence on a new line as it 
  decreases maintenance and makes diffs more readable.

.. title:: Information for uap Developers

..
  This document describe different issues concerning the development of **uap**.


*************
Documentation
*************

The official documentation is hosted on readthedocs (add link). But as
developer you need to be able to create the documentation locally. So let's
focus at first on:

Creating the Documentation Locally
==================================



Prerequisites
-------------

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

The documentation can be build by::

  $ cd doc
  $ make html

This should do the trick.
