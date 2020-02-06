#!/bin/bash

set -e

##################
# Submit Options #
##################

#SBATCH #{SUBMIT_OPTIONS}          # Arbitrary submit options

############
# Commands #
############

array_jobs=(#{ARRAY_JOBS})

#{PRE_JOB_COMMAND}

#{COMMAND}

#{POST_JOB_COMMAND}
