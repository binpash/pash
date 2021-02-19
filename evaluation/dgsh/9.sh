#!/bin/bash
#
# SYNOPSIS C/C++ symbols that should be static
# DESCRIPTION
# Given as an argument a directory containing object files, show which
# symbols are declared with global visibility, but should have been
# declared with file-local (static) visibility instead.
# Demonstrates the use of dgsh-capable comm (1) to combine data from
# two sources.

# Find object files
mkfifo a b c d

# List all defined (exported) symbols
cat a | awk 'NF == 3 && $2 ~ /[A-Z]/ {print $3}' | sort > c&

# List all undefined (imported) symbols
cat b | awk '$1 == "U" {print $2}' | sort  > d &

find "$1" -name \*.o | xargs nm | tee a b  > /dev/null 
comm -23 c d
#comm -23
rm a b c d
