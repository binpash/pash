#!/bin/bash

# $input="${1}"
# $output="${2}"
# $temp="${3}"

touch "$3"

cat "$1" > "$3"
cat "$3" > "$2"
