#!./python_env/bin/python

import sys
sys.path.append('./include')
import pipeline
import yaml

run_mode = pipeline.Pipeline.run_modes.FULL
if len(sys.argv) > 1:
    if sys.argv[1] == '--dry-run':
        run_mode = pipeline.Pipeline.run_modes.DRY_RUN
    elif sys.argv[1] == '--test-run':
        run_mode = pipeline.Pipeline.run_modes.TEST_RUN

p = pipeline.Pipeline(run_mode)

p.print_tasks()
