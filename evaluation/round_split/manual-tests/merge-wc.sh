#!/bin/bash

A="paste -d '+' "
for i in "$@"; do
  # cat "$i" | tr -s ' '  '\n' | tail -n +2
  A="$A <(cat $i | tr -s ' '  '\n' | tail -n +2) "
done
A="$A | head -n +3 | bc | tr -s '\n'  ' ' |  sed 's/$/\ /'"

eval $A

