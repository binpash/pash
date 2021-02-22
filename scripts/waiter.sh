#!/usr/bin/env bash

## This script runs a command in the background
## and waits until the command closes its stdout.
## When it does, it writes to a designated file to
## show command completion.

## If you want to run the waiter on the background (which is the main reason to use it),
## run it using `nohup waiter.sh ... &`.

if [ $# -lt 3 ]; then
 echo "Not enough arguments!"
 exit 1
fi

completion_file=$1
output_log=$2
command="${@:3}"

rm -f "$completion_file"

pipe_name="$(mktemp -u)"
mkfifo "$pipe_name"
$command > "$pipe_name" &
cat $pipe_name > "$output_log"

echo "done" > "$completion_file"
