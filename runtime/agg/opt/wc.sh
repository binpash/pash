#!/usr/bin/env bash

A="paste -d '+' "
for i in "$@"; do
  # cat "$i" | tr -s ' '  '\n' | tail -n +2
  A="$A <(cat $i | tr -s ' '  '\n' | tail -n +2) "
done
A="$A | bc | tr -s '\n'  ' ' | sed 's/^/   /' | sed 's/$/\n /'"

eval "$A"
