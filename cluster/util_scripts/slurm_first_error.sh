#! /usr/bin/env bash

# Pass space seperated slurm job ids to this script to recive
# the stderr produced by the first job that failed.

shopt -s lastpipe

for jobid in "$@"; do
    scontrol show job -o $jobid
done | grep "JobState=FAILED .* StdErr=" |
    sed -E 's/^JobId=([0-9]*) .*EndTime=([^ ]*) .*/\2 \1/g' |
    sort -u | head -n1 | sed 's/[^ ]* //' |
    read stderrJob


if [[ -z "$stderrJob" ]]; then
    >&2 printf 'No failed job could be found with scontrol for the jobs: %s\n' "$*"
else
    printf 'The first job that failed was %s.\n' "$stderrJob"
    cmd="scontrol show job $stderrJob"
    printf '\e[31mOutput of: %s\e[0m\n' "$cmd"
    printf -v esc '\e'
    eval "$cmd" | sed "s/\([^ ,]*\)=/\x0${esc}[34m\1: ${esc}[0m/g"
    stderrFile="$(eval "$cmd -o" | sed -E 's/.* StdErr=(.*) StdIn=.*/\1/')"
    printf '\e[31mContent of: %s\e[0m\n' "$stderrFile"
    cat "$stderrFile"
fi
