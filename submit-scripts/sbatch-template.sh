#!/bin/bash
##SBATCH --job-name=md5sum
##SBATCH --output="my_test_job_%j.txt"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=#{CORES}


#{COMMAND}




