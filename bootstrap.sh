#!/usr/bin/env bash

# This script creates a virtual Python environment in ./python-env with all required libraries set up.
# To use it, call Python at ./python-env/bin/python.

virtualenv python_env
./python_env/bin/pip install pyyaml numpy biopython psutil
gcc -o tools/cat4m tools/cat4m.c
git submodule update --init
cd tools/pigz4m; make; cd ../..

