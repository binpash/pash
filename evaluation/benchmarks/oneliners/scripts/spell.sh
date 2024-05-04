#!/bin/bash
# Calculate mispelled words in an input
# https://dl.acm.org/doi/10.1145/3532.315102

cd "$(dirname "$0")" || exit 1

[ -z "$PASH_TOP" ] && {
  echo "PASH_TOP not set, maybe $(git rev-parse --show-toplevel)?"
  exit
}

IN=${IN:-$PASH_TOP/evaluation/benchmarks/oneliners/inputs/1G.txt}
DICT=${DICT:-$PASH_TOP/evaluation/benchmarks/oneliners/inputs/dict.txt}

cat "$IN" |
    iconv -f utf-8 -t ascii//translit | # remove non utf8 characters
    col -bx |                           # remove backspaces / linefeeds
    tr -cs A-Za-z '\n' |
    tr A-Z a-z |                        # map upper to lower case
    tr -d '[:punct:]' |                 # remove punctuation
    sort |                              # put words in alphabetical order
    uniq |                              # remove duplicate words
    comm -23 - "$DICT"                  # report words not in dictionary
