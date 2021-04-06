#!/bin/bash
# tag: highlight misspelled words
# from: https://www2.dmst.aueb.gr/dds/sw/dgsh/#spell-highlight
set -e
export LC_ALL=C
IN=${VOC:-/usr/share/dict/words}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/dgsh/input}
sort /usr/share/dict/words | cat > ${OUT}/sort_dict
mkfifo ${OUT}/a ${OUT}/sort_IN
cat ${OUT}/a |	tr -cs A-Za-z \\n | tr A-Z a-z | sort -u | cat > ${OUT}/sort_IN &
# tee process
cat ${IN} | tee ${OUT}/a  > /dev/null &
comm -23 ${OUT}/sort_IN ${OUT}/sort_dict | cat > ${OUT}/bad && grep --fixed-strings --file=- --ignore-case --color --word-regex --context=2 < ${OUT}/bad ${IN}

rm -f ${OUT}/a ${OUT}/sort_IN ${OUT}/sort_dict ${OUT}/bad
