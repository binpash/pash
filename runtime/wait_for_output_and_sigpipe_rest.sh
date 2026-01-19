#!/usr/bin/env bash

## TODO: Give it the output pid as an argument
wait "$@"

## TODO: This only works if there is only a single output node 
##         (and a single node given as an argument to `wait`).
export internal_exec_status=$?

# It is assumed that $distro is set when this is called.

# Note: We need the || true after the grep so that it doesn't exit with error if it finds nothing.


# This value may contains multiple pids as a whitespace-separated string, and
# we must split it as multiple pids by shell's field splitting.
# shellcheck disable=SC2086
(unset IFS; > /dev/null 2>&1 kill -SIGPIPE $pids_to_kill || true)

##
## Old way of waiting, very inefficient.
##

# now do different things depending on distro

## TODO: Delete this since it is very costly
# case "$distro" in
#     freebsd*)  
#         # not sure at all about this one
#         pids_to_kill="$(ps -efl $BASHPID |awk '{print $1}' | { grep -E '[0-9]' || true; } )"
#         ;;
#     *)
#         pids_to_kill="$(ps --ppid $BASHPID |awk '{print $1}' | { grep -E '[0-9]' || true; } )"
#         ;;
# esac
# pids_to_kill=""

## TODO: Maybe send a signal to all pids at once
# for pid in $pids_to_kill
# do
#     # wait $pid
#     (> /dev/null 2>&1 kill -SIGPIPE $pid || true)
# done
