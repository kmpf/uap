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

virtualenv python_env
./python_env/bin/pip install pyyaml psutil biopython 
./python_env/bin/easy_install -f http://biopython.org/DIST/ biopython

# Installation of NumPy on EVE #

# Backup LDFLAGS and CPPFLAGS values 
export LDFLAGS_BAK=$LDFLAGS && export CPPFLAGS_BAK=$CPPFLAGS
# unset environment variables
unset LDFLAGS && unset CPPFLAGS
./python_env/bin/pip install numpy
export LDFLAGS=$LDFLAGS_BAK && export CPPFLAGS=$CPPFLAGS_BAK
unset LDFLAGS_BAK && unset CPPFLAGS_BAK

module load freetype/2.5.5-1
./python_env/bin/pip install matplotlib

git submodule update --init --recursive

