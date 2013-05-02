#!/bin/bash

#$ -pe smp #{CORES}           # request the "smp" parallel environment and up to #{CORES} slots
#$ -cwd                       # Execute the job from the current working directory. This switch will activate Sun Grid Engine’s path aliasing facility, if the corresponding configuration files are present (see sge_aliases(5)).
#$ -V                         # Specifies that all environment variables active within the qsub utility be exported to the context of the job.
#$ -S /bin/bash               # shell to be used
#$ -m as                     # e-mail notification at following events: End of job, Abortion of job and Suspension of job
#$ -M #{EMAIL}                # e-mail notification
#$ -l h_rt=96:00:00           # runtime of your job
#$ -l m_core=#{CORES}         # use only nodes where all #{CORES} cores are available

# loading important environment stuff
[ -r ~/.bashrc ]       && . ~/.bashrc
[ -r ~/.bash_profile ] && . ~/.bash_profile

##
## does the actual work
##
source /etc/profile.d/000-modules.sh    # makes the module system available to your job
module load python/2.7.3-1-virtual
module load pigz

#{COMMAND}
