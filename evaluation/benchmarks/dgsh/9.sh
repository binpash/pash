#!/bin/bash
#
# SYNOPSIS C/C++ symbols that should be static
# DESCRIPTION
# Given as an argument a directory containing object files, show which
# symbols are declared with global visibility, but should have been
# declared with file-local (static) visibility instead.
# Demonstrates the use of dgsh-capable comm (1) to combine data from
# two sources.

# tag: symbols that should be static
set -e
IN=${IN:-/usr/local/bin}
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

find ${IN} -type f | xargs nm | tee a b  > /dev/null 
comm -23 c d
rm a b c d
