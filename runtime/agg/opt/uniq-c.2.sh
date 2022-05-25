#!/usr/bin/env bash

END_OF_1=$(tail -n 1 "$1")
END_NUM=$(echo "$END_OF_1" | grep -E -o '^[ ]*[0-9]*[ ]*' | tr -d "[:space:]")
END_WORD=$(echo "$END_OF_1" | sed 's/^[ ]*[0-9]*[ ]*//g')

START_OF_2=$(head -n 1 "$2")
START_NUM=$(echo "$START_OF_2" | grep -E -o '^[ ]*[0-9]*[ ]*' | tr -d "[:space:]")
START_WORD=$(echo "$START_OF_2" | sed 's/^[ ]*[0-9]*[ ]*//g')

if [[ $START_WORD == "$END_WORD" ]]; then
  TOTAL_NUM=$((START_NUM + END_NUM))
  sed '$d' "$1"
  printf "%7s %s\n" "$TOTAL_NUM" "$START_WORD"
  sed '1d' "$2"
else
  cat "$1" "$2"
fi
