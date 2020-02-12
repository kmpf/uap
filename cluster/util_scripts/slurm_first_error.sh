#! /usr/bin/env bash

# Pass space seperated slurm job ids to this script to recive
# the stderr produced by the first job that failed.

shopt -s lastpipe

for jobid in "$@"; do
    scontrol show job -o $jobid
done | grep -v "JobState=COMPLETED" |
    grep -v "JobState=RUNNING" |
    grep -v "JobState=PENDING" |
    sed -E 's/^JobId=([0-9]*) .*EndTime=([^ ]*) .*/\2 \1/g' |
    sort -u | head -n1 | sed 's/[^ ]* //' |
    read stderrJob


if [[ -z "$stderrJob" ]]; then
    >&2 printf 'All found jobs were "PENDING", "RUNNING" or "COMPLETED".\n'
else
    printf 'The first job that is not "PENDING", "RUNNING" or "COMPLETED" is %s.\n' "$stderrJob"
    cmd="scontrol show job $stderrJob"
    printf '\e[31mOutput of: %s\e[0m\n' "$cmd"
    printf -v esc '\e'
    eval "$cmd" |
        grep -B99 -m1 "^$" | # in case JobId==ArrayJobId
        sed "s/\([^ ,]*\)=/\x0${esc}[34m\1: ${esc}[0m/g"
    stderrFile="$(eval "$cmd -o" | head -n1 |
        sed -E 's/.* StdErr=(.*) StdIn=.*/\1/')"
    printf '\e[31mLast 30 lines of: %s\e[0m\n' "$stderrFile"
    tail -n30 "$stderrFile"
fi
