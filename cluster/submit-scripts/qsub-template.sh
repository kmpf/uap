#!/bin/bash

set -e

##################
# Submit Options #
##################

#$ #{SUBMIT_OPTIONS}          # User defined submit options

############
# Commands #
############

array_jobs=(#{ARRAY_JOBS})

#{PRE_JOB_COMMAND}

#{COMMAND}

#{POST_JOB_COMMAND}
