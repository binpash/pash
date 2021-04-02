#!/usr/bin/env dgsh
# not working
export LC_ALL=C
INPUT=${INPUT:-/usr/share/dict/words}
OUTPUT=${OUTPUT:-$PASH_TOP/evaluation/dgsh/output}
sort /usr/share/dict/words | cat > ${OUTPUT}/sort_dict
mkfifo ${OUTPUT}/a ${OUTPUT}/sort_input
cat ${OUTPUT}/a |	tr -cs A-Za-z \\n | tr A-Z a-z | sort -u | cat > ${OUTPUT}/sort_input &
# tee process
cat ${INPUT} | tee ${OUTPUT}/a  > /dev/null &
comm -23 ${OUTPUT}/sort_input ${OUTPUT}/sort_dict | cat > ${OUTPUT}/bad && grep --fixed-strings --file=- --ignore-case --color --word-regex --context=2 < ${OUTPUT}/bad ${INPUT}

rm -f ${OUTPUT}/a ${OUTPUT}/sort_input ${OUTPUT}/sort_dict ${OUTPUT}/bad
