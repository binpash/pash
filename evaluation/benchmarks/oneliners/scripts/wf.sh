#!/bin/bash
# Calculate the frequency of each word in the document, and sort by frequency

cd "$(dirname "$0")" || exit 1

[ -z "$PASH_TOP" ] && {
  echo "PASH_TOP not set, maybe $(git rev-parse --show-toplevel)?"
  exit
}

IN=${IN:-$PASH_TOP/evaluation/benchmarks/oneliners/inputs/1G.txt}

cat "$IN" | tr -cs A-Za-z '\n' | tr A-Z a-z | sort | uniq -c | sort -rn
