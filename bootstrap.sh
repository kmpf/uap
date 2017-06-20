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
./python_env/bin/pip install pyyaml numpy biopython psutil

# create a virtualenv for a different python module
# in this case its essential for using the multiprocessing functions
# (e.g. in tool 'segemehl_2017_reformatCigar.py'

virtualenv -p /usr/local/python/2.7.13-1/bin/python python_2_7.13-1_env
