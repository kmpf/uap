#!/bin/bash

##$ -pe smp #{CORES}           # request the "smp" parallel environment and up
#                              # to #{CORES} slots
##$ -cwd                       # Execute the job from the current working
#                              # directory. This switch will activate Sun Grid
#                              # Engine’s path aliasing facility, if the
#                              # corresponding configuration files are present
#                              # (see sge_aliases(5)).
##$ -S /bin/bash               # shell to be used
##$ -m as                      # E-mail notification on aborting and suspending
#                              # jobs
##$ -M #{EMAIL}                # E-mail notification address
##$ -l h_rt=#{TIME}            # runtime of your job
##$ -l h_vmem=#{MEMORY}        # memory per core for your job

#$ #{SUBMIT_OPTIONS}          # User defined submit options

## Execute command

##{MODULE_LOAD}

#{PRE_JOB_COMMAND}

#{COMMAND}

#{POST_JOB_COMMAND}

##{MODULE_UNLOAD}
