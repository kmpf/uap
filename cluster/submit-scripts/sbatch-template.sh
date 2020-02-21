#!/bin/bash

##################
# Submit Options #
##################

#SBATCH #{SUBMIT_OPTIONS}          # Arbitrary submit options

############
# Commands #
############

set -e

config=$(mktemp)
cat << 'EOF' > $config

#{UAP_CONFIG}

EOF
exec 123< $config
rm $config

array_jobs=(#{ARRAY_JOBS})

#{PRE_JOB_COMMAND}

#{COMMAND}

#{POST_JOB_COMMAND}
