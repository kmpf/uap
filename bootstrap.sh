#!/usr/bin/env bash

# This script creates a virtual Python environment in ./python-env with all required libraries set up.
# To use it, call Python at ./python-env/bin/python.

virtualenv python_env
./python_env/bin/pip install textile pyyaml
