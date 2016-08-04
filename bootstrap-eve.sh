#!/usr/bin/env bash

# This script creates a virtual Python environment in ./python-env with all required libraries set up.
# To use it, call Python at ./python-env/bin/python.

# Test if everything is available
virtualenv --version >/dev/null 2>&1 || missing=$missing"[bootstrap.sh] Python virtualenv is required but it's not installed.\n"
gcc --version >/dev/null 2>&1 || missing=$missing"[bootstrap.sh] GCC is required but it's not installed.\n"
git --version >/dev/null 2>&1 || missing=$missing"[bootstrap.sh] GIT is required but it's not installed.\n"

# Report errors if any
if [[ -n "$missing" ]]; then
    echo -e ${missing%\\n}
    echo "[bootstrap.sh] Please install missing programs."
    exit 1
fi

################################
# Installation of NumPy on EVE #
################################

# Backup LDFLAGS and CPPFLAGS values 
export LDFLAGS_BAK=$LDFLAGS && export CPPFLAGS_BAK=$CPPFLAGS
# unset environment variables
unset LDFLAGS && unset CPPFLAGS
# Install numpy
./python_env/bin/pip install numpy
# Reset LDFLAGS and CPPFLAGS values
export LDFLAGS=$LDFLAGS_BAK && export CPPFLAGS=$CPPFLAGS_BAK
# Remove backup variables
unset LDFLAGS_BAK && unset CPPFLAGS_BAK

##############################################
# Installation of PyYAML, psutil, BioPython  #
##############################################

# Create virtual python environment named 'python_env'
virtualenv python_env
# Install required libraries (pyyaml,psutil,biopython) into the virtual python environemnt
./python_env/bin/pip install pyyaml psutil biopython 
./python_env/bin/easy_install -f http://biopython.org/DIST/ biopython

##############################
# Installation of matplotlib #
##############################

# matplotlib requires freetype
module load freetype/2.5.5-1
./python_env/bin/pip install matplotlib
