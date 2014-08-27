#!/bin/bash

#$ -pe smp #{CORES}           # request the "smp" parallel environment and up to #{CORES} slots
#$ -cwd                       # Execute the job from the current working directory. This switch will activate Sun Grid Engine’s path aliasing facility, if the corresponding configuration files are present (see sge_aliases(5)).
##$ -V                         # Specifies that all environment variables active within the qsub utility be exported to the context of the job.
#$ -S /bin/bash               # shell to be used
#$ -m as                      # e-mail notification on aborting and suspending jobs
#$ -M #{EMAIL}                # e-mail notification address
#$ -l h_rt=96:00:00           # runtime of your job
#$ -l h_vmem=5G               # memory per core for your job
#$ -l m_core=#{CORES}         # use only nodes where all #{CORES} cores are available

# loading important environment stuff
[ -r ~/.bashrc ]       && . ~/.bashrc
[ -r ~/.bash_profile ] && . ~/.bash_profile

##
## does the actual work
##
source /etc/profile.d/000-modules.sh    # makes the module system available to your job
module load emacs/24.3-2_gcc_4.8.1_CentOS6
    
# Load language modules and support modules
module load R/3.1.0-1_gcc_4.8.1_CentOS6
module load perl/5.20/1-1_CentOS6
module load python/2.7.6-2_gcc_4.8.1_CentOS6-virtual

# Load bioinformatics modules
module load bamtools/2.3.0-1_CentOS6
module load bedtools/2.20.1-1_CentOS6
module load samtools/0.1.19-1_CentOS6
module load tophat/2.0.12-1_CentOS6
module load picard-tools/1.117-1_CentOS6
module load casava/1.8.2-4_CentOS6
module load fastqc/0.11.2-1_CentOS6
module load cutadapt/1.4.1-1_CentOS6
module load segemehl/0.1.9-1_CentOS6
    
# Load support modules
module load git/1.9.2-2_gcc_4.8.1_CentOS6
module load lapack/3.5.0-1_gcc_4.8.1_CentOS6


#{COMMAND}
