<!--
  This is the documentation for uap. Please keep lines under
  80 characters if you can and start each sentence on a new line as it 
  decreases maintenance and makes diffs more readable.

  See https://github.github.com/gfm/ for GitHub flavored Markdown
-->

# uap -- Robust, Consistent, and Reproducible Data Analysis

**Authors:**

Christoph Kämpf, Michael Specht, Alexander Scholz, Sven-Holger Puppel,
Gero Doose, Kristin Reiche, Jana Schor, Jörg Hackermüller

**Description:**

  **uap** executes, controls and keeps track of the analysis of large data sets.
  It enables users to perform robust, consistent, and reproducible data analysis.
  **uap** encapsulates the usage of (bioinformatic) tools and handles data flow
  and processing during an analysis.
  Users can use predefined or self-made analysis steps to create custom analysis.
  Analysis steps encapsulate best practice usages for bioinformatic software
  tools.
  **uap** focuses on the analysis of high-throughput sequencing (HTS) data.
  But its plugin architecture allows users to add functionality, such that
  it can be used for any kind of large data analysis.

**Usage:**

  **uap** is a command-line tool, implemented in Python.
  It requires a user-defined configuration file, which describes the analysis,
  as input.
  
**Supported Platforms:**

  * Unix-like operating systems.
  * High Performance Compute (HPC) cluster systems such as
    [UGE](http://www.univa.com/products/), [OGE/SGE](http://www.univa.com/oracle)
    and [SLURM](http://slurm.schedmd.com/).
  * see [Tested Platforms](./platforms.md) for detailed information
    
**Important Information**

  **uap** does **NOT** include all tools necessary for the data analysis.
  It expects that the required tools are **already installed**.

## Table of contents

[Consistency](./introduction.md#consistency)

1. [Introducing uap](./introduction.md)
2. [Quick Start uap](./how-to.md)
3. [](./installation.md)
4. [](./configuration.md)
5. [](./interaction.md)
6. [](./extension.md)
7. [](./platforms.md)
8. [](./annotation.md)
9. [](./steps.md)
10. [](./api.md)
11. [](./developers.md)

## Remarks

This documentation uses [GitHub Pages]() and is written in [Markdown]().
