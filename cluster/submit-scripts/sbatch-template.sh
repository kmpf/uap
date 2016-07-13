#!/bin/bash

#SBATCH --cpus-per-task=#{CORES}   # Request #{CORES} for this task
#SBATCH --ntasks=1
#SBATCH #{SUBMIT_OPTIONS}          # Arbitrary submit options

## Execute command

#{COMMAND}




