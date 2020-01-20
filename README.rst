uap -- Universal Analysis Pipeline
==================================

The **uap** package is a framework to configure, run, and control
large data multi-step analyses.
Its main focus is on the analysis of high-throughput sequencing data.

The aim of this data processing pipeline is to enable robust and straightforward
bioinformatics data evaluation.
It is implemented in Python, runs under GNU/Linux and can be controlled from the
command-line interface.
Although the primary focus is the evaluation of sequencing data, its design
allows for a variety of other applications.


Documentation
=============

The documentation of **uap** is available as `Giltab Page <https://onebutton.ribogitpages.izi.fraunhofer.de/uap/>`_.

Testing
=======

In order to use the testing repo [uap_test](https://ribogit.izi.fraunhofer.de/oneButton/uap_test)
you hav to change the file .gitmodules
by replacing `url = ../uap_test.git` with `url = git@ribogit.izi.fraunhofer.de:oneButton/uap_test.git`.
The entry for `uap_test` looks like this:

```
[submodule "uap_test"]
        path = uap_test
        url = git@ribogit.izi.fraunhofer.de:oneButton/uap_test.git
```

Then you can run `git submodule sync && git checkout HEAD -- .gitmodules && git submodule update`
to have the testing repo in `uap_test`. Please consult tutorials for submodules for further info.
E.g., https://git-scm.com/book/en/v2/Git-Tools-Submodules.
