#!/bin/sh

# FIXME: the following expansion would not have worked
# (https://github.com/andromeda/pash/issues/218)
# cd "$(dirname $0)"

[ -z $PASH_TOP ] && { 
  echo "PASH_TOP not set, maybe $(git rev-parse --show-toplevel)?"
  exit
}
FILE="$PASH_TOP/evaluation/intro/input/100M.txt"
DICT="$PASH_TOP/evaluation/intro/input/sorted_words"

cat "$FILE" | tr A-Z a-z | tr -cs A-Za-z '\n' | sort | uniq | comm -13 $DICT -
