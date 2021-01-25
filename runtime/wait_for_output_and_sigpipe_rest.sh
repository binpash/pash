#!/bin/bash

## TODO: Give it the output pid as an argument

wait $@
pids_to_kill="$(ps --ppid $$ |awk '{print $1}' | grep -E '[0-9]')"
for pid in $pids_to_kill
do
        (> /dev/null 2>&1 kill -SIGPIPE $pid || true)
done
