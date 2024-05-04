#!/bin/bash
# Sort input

cd "$(dirname "$0")" || exit 1

SIZE=500M

IN=input/$SIZE.txt

cat $IN | sort
