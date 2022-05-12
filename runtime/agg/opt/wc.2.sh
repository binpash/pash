#!/usr/bin/env bash

# Part of a distributed-`wc` wrapper, merging two `wc` results
# FIXME needs correct padding

paste -d '+'
    <(cat "$1" |
      wc |
      tr -s ' '  '\n' |
      tail -n +2)
    <(cat "$2" |
      wc |
      tr -s ' '  '\n' |
      tail -n +2) |
  bc |
  tr -s '\n'  ' ' |
  sed 's/^/   /' |
  sed 's/$/\ /'
