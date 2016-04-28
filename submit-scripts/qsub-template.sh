#!/bin/bash

#$ -pe smp #{CORES}           # request the "smp" parallel environment and up to #{CORES} slots
#$ -cwd                       # Execute the job from the current working directory. This switch will activate Sun Grid Engine’s path aliasing facility, if the corresponding configuration files are present (see sge_aliases(5)).
##$ -V                         # Specifies that all environment variables active within the qsub utility be exported to the context of the job.
#$ -S /bin/bash               # shell to be used
#$ -m as                      # e-mail notification on aborting and suspending jobs
#$ -M #{EMAIL}                # e-mail notification address
#$ -l h_rt=#{TIME}            # runtime of your job
#$ -l h_vmem=#{MEMORY}        # memory per core for your job

## Execute command

#{COMMAND}
