# count how many times each file type exist in a directory
INPUT=${INPUT:-$PASH_TOP/evaluation/aliases/input/}
OUTPUT=${OUTPUT:-$PASH_TOP/evaluation/aliases/output}
find $INPUT -type f | while read f; do echo ""${f##*.}""; done | sed ""/^\s*$/d"" | sort | uniq -c | sort -rn  > $OUTPUT/get_type_count_res
