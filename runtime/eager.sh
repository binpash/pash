#!/usr/bin/env bash

input=${1?"ERROR: Eager: No input file given"}
output=${2?"ERROR: Eager: No output file given"}
intermediate_file=${3?"ERROR: Eager: No intermediate file given"}

# Set a default DISH_TOP in this directory if it doesn't exist
PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}

cleanup()
{
    kill -SIGPIPE $eager_pid > /dev/null 2>&1
}
trap cleanup EXIT

$PASH_TOP/runtime/eager "$input" "$output" "$intermediate_file" &
eager_pid=$!
wait $eager_pid
rm "$intermediate_file"
