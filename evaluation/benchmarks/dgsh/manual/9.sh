#!/bin/bash
# from: https://www2.dmst.aueb.gr/dds/sw/dgsh/#static-functions
# tag: symbols that should be static
set -e
IN=${BIN:-/usr/local/bin}
#OUTPUT=${OUTPUT:-$PASH_TOP/evaluation/benchmarks/dgsh/output}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/dgsh/input}
cd ${OUT}
# Find object files
rm -f a b c d
mkfifo a b c d

# List all defined (exported) symbols
cat a | awk 'NF == 3 && $2 ~ /[A-Z]/ {print $3}' | sort > c&

# List all undefined (imported) symbols
cat b | awk '$1 == "U" {print $2}' | sort  > d &

find ${IN} -type f -exec grep -IL . "{}" \; | xargs nm | tee a b  > /dev/null 
comm -23 c d
rm a b c d
