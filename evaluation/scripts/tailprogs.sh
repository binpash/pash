#!/bin/bash

# A bash script for finding the 10 longest scripts 
# (TODO: `group_by` script type?)

# From "Wicked Cool Shell Scripts", 2nd Ed., pg. 7

# Data:
# Assumes a full list of commands
#
# # simple, from a single dir:
# echo "$(
#   ls /usr/bin/*
# )" > all_cmds.txt
# 
# # Or more complicated, from $PATH:
#
# echo "$(
#   case "$PATH" in
#     (*[!:]:) PATH="$PATH:" ;;
#   esac
#   
#   set -f; IFS=:
#   for dir in $PATH; do
#     set +f
#     [ -z "$dir" ] && dir="."
#     for file in "$dir"/*; do
#       if [ -x "$file" ] && ! [ -d "$file" ]; then
#         printf '%s = %s\n' "${file##*/}" "$file"
#       fi
#     done
#   done
# )" > ./input/allcmds.txt

IN=./input/cmds10x.txt
OUT=./output/out.txt

ls /usr/bin/* > $IN

cat $IN |
  xargs file |
  grep "shell script" |
  cut -d: -f1 |
  xargs wc -l |
  sort -rn |
  head -n 25 > $OUT
