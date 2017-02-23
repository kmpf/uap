#!/bin/bash

##################
# Submit Options #
##################

#SBATCH #{SUBMIT_OPTIONS}          # Arbitrary submit options

############
# Commands #
############

#{PRE_JOB_COMMAND}

#{COMMAND}

#{POST_JOB_COMMAND}
