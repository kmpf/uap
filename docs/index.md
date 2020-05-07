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

## Table of Content

1. [Introducing uap](./introduction.md)
2. [Quick Start](./how-to.md)
3. [Installation](./installation.md)
4. [Analysis Configuration File](./configuration.md)
5. [Command-Line Usage](./interaction.md)
6. [Add New Functionality](./extension.md)
7. [Tested Platforms](./platforms.md)
8. [Annotation Files](./annotation.md)
9. [Available Steps](./steps.md)
10. [API Documentation](./api.md)
11. [Information for developers](./developers.md)

<dl>
  <dt>First Term</dt>
  <dd>This is the definition of the first term.</dd>
  <dt>Second Term</dt>
  <dd>This is one definition of the second term. </dd>
  <dd>This is another definition of the second term.</dd>
</dl>
