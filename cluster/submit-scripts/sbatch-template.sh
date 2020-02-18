#!/bin/bash

set -e

##################
# Submit Options #
##################

#SBATCH #{SUBMIT_OPTIONS}          # Arbitrary submit options

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
