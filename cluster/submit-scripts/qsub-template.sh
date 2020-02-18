#!/bin/bash

set -e

##################
# Submit Options #
##################

#$ #{SUBMIT_OPTIONS}          # User defined submit options

############
# Commands #
############

config=$(mktemp)
cat << EOF > $config

#{UAP_CONFIG}

EOF

array_jobs=(#{ARRAY_JOBS})

#{PRE_JOB_COMMAND}

#{COMMAND}

#{POST_JOB_COMMAND}
