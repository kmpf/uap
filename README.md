# uap -- Universal Analysis Pipeline

## Authors

Christoph Kämpf, Michael Specht, Alexander Scholz, Sven-Holger Puppel, Gero Doose, Kristin Reiche, Jana Schor, Jörg Hackermüller

[uap: reproducible and robust HTS data analysis. BMC Bioinformatics 20, 664 (2019)](https://doi.org/10.1186/s12859-019-3219-1)

## Introduction

The **uap** package is a framework to configure, run, and control
large data multi-step analyses.
Its main focus is on the analysis of high-throughput sequencing data.

The aim of this data processing pipeline is to enable robust and straightforward
bioinformatics data evaluation.
It is implemented in Python, runs under GNU/Linux and can be controlled from the
command-line interface.
Although the primary focus is the evaluation of sequencing data, its design
allows for a variety of other applications.

## About this Repository

This repository contains the development status of **uap** at Fraunhofer IZI.
It is based on the **uap** repository as published in [Kämpf, C., Specht, M.,Scholz, A. et al. uap: reproducible and robust HTS data analysis. BMC Bioinformatics 20, 664 (2019)](https://doi.org/10.1186/s12859-019-3219-1), which is located [here](https://github.com/yigbt/uap).
 
This version [v2.0.0rc2](https://github.com/fraunhofer-izi/uap/releases/tag/v2.0.0rc2) contains the following changes (for a complete list see the [CHANGELOG](CHANGELOG.md)):

* code conversion from Python2 to Python3
* improved user interaction
* enhanced error detection for configuration
* validation of existing results by using annotation as configuration and recalculation of SHA256
* `status --details` completely lists errors or changes caused by adaptation of the configuration 
* enhanced detection of changes (software version, output files, sha256 of results (optional))
* improved error-management
* removed checksum suffix in output directories
* extended backward-compatible connection-management
* Source_controller step to check input data
* no need to configure `uap` internal scripts, GNU coreutils and `lmod`
* improved job-execution (signal handling, array jobs, enhanced logging, changes to configuration do not impact running jobs)
* processes are executed in temporary directories
* error fixing and code improvement

Please note, the version [v2.0.0rc2](https://github.com/fraunhofer-izi/uap/releases/tag/v2.0.0rc2) of **uap** requires Python >= 3.5 and is only tested on [SLURM](https://slurm.schedmd.com/documentation.html).

## Contacts

* Christoph Kämpf, [christoph.kaempf@izi.fraunhofer.de](mailto:christoph.kaempf@izi.fraunhofer.de)
* Kristin Reiche, [kristin.reiche@izi.fraunhofer.de](mailto:kristin.reiche@izi.fraunhofer.de)
* Jana Schor, [jana.schor@ufz.de](mailto:jana.schor@ufz.de)
* Jörg Hackermüller, [joerg.hackermueller@ufz.de](mailto:joerg.hackermueller@ufz.de)

Helmholtz Centre for Environmental Research - UFZ</br>
Permoserstr. 15, 04318 Leipzig, Germany

Fraunhofer Institute for Cell Therapy and Immunology (IZI)</br>
Perlickstraße 1, 04103 Leipzig, Germany

## Main contributors:

* Christoph Kämpf
* Dominik Otto
* Michael Specht
* Alexander Scholz
* Sven-Holger Puppel
* Gero Doose
* Kristin Reiche
* Jana Schor
* Jörg Hackermüller

## Copyright Notice

Copyright (C) 2011  - 2020, *Helmholtz-Zentrum für Umweltforschung GmbH – UFZ*
and *Fraunhofer-Gesellschaft zur Förderung der angewandten Forschung e.V.*
acting on behalf of its *Fraunhofer Institute for Cell Therapy and Immunology
(IZI)*. All rights reserved.

### The code is a property of:

*Helmholtz-Zentrum für Umweltforschung GmbH – UFZ*</br>
Registered Office: Leipzig</br>
Registration Office: Amtsgericht Leipzig</br> 
Trade Register Nr. B 4703</br>
Chairman of the Supervisory Board: MinDirig'in Oda Keppler</br>
Scientific Director: Prof. Dr. Georg Teutsch</br>
Administrative Director: Dr. Sabine König

and 

*Fraunhofer-Gesellschaft zur Foerderung der angewandten Forschung e.V.*</br>
Hansastraße 27c, 80686 München</br>
acting on behalf of its</br>
*Fraunhofer Institute for Cell Therapy and Immunology (IZI)*</br>
Perlickstraße 1, 04103 Leipzig

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program; if not, it can be found at the end of this document or see at
<https://www.gnu.org/licenses/gpl-3.0.de.html>.

## Redistribution

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, the
  list of conditions for redistribution and modification as well as the
  following GNU General Public License.
* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions, the following GNU General Public License and the
  modification conditions in the documentation and/or other materials provided
  with the distribution.
* Neither the name of *Helmholtz-Zentrum für Umweltforschung GmbH – UFZ*, the name
  of the *Fraunhofer Institute for Cell Therapy and Immunology (IZI)* nor the
  names of its contributors may be used to endorse or promote products derived
  from this software without specific prior written permission.

## Modification

If software is modified to produce derivative works, such modified software
should be clearly marked, so as not to confuse it with the version available
from *Helmholtz-Zentrum für Umweltforschung GmbH – UFZ* or *Fraunhofer
Institute for Cell Therapy and Immunology (IZI)*.
