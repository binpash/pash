#!/bin/bash

# This classification  is wrt to their operation, not  its input---i.e., whether
# input contains the use of fs identifiers (identifiers can be fs or not fs)

# We need to think about how to translate DFS commands
# What is a distributed fs? Directories are simply keys?

# everything else (i.e., side-effectful) just needs to be converted to location independent commands

p="../c_stats/"
A=${1:-${p}posix.txt}
B=${2:-${p}coreutils.txt}

# Take commands that are shared and use existing distributability descriptions
comm -12 <(cat $A | grep 'Mandatory' | cut -d ' ' -f 1 | sort ) <( cut -d ' ' -f 1 $B | sort) |
  sed s/^/\^/ |
  xargs -n 1 -I {} grep -w {} ./coreutils.txt |
  sort -b -k2,2 -k1,1 # > posix_mandatory1.txt # commenting out this redirection will overwrite!

# Analyze mandatory commands not in the second, and not built-ins
comm -23 <(cat $A | grep 'Mandatory' | cut -d ' ' -f 1 | sort ) <( cut -d ' ' -f 1 $B | sort) | 
  comm -23 - <(cat ../c_stats/builtins.txt | sed 's/  */ /g' | cut -d ' ' -f 1 | sort) |
  sed s/^/\^/ |
  xargs -n 1 -I {} grep -w {} $A |
  sed s/Mandatory// |
  sort -b -k2,2 -k1,1 # > posix_mandatory2.txt # commenting out this redirection will overwrite!
