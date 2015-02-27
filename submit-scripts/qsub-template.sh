#!/bin/bash

#$ -pe smp #{CORES}           # request the "smp" parallel environment and up to #{CORES} slots
#$ -cwd                       # Execute the job from the current working directory. This switch will activate Sun Grid Engine’s path aliasing facility, if the corresponding configuration files are present (see sge_aliases(5)).
##$ -V                         # Specifies that all environment variables active within the qsub utility be exported to the context of the job.
#$ -S /bin/bash               # shell to be used
#$ -m as                      # e-mail notification on aborting and suspending jobs
#$ -M #{EMAIL}                # e-mail notification address
#$ -l h_rt=96:00:00           # runtime of your job
#$ -l h_vmem=8G               # memory per core for your job
#$ -l m_core=#{CORES}         # use only nodes where all #{CORES} cores are available

## Execute command

#{COMMAND}
