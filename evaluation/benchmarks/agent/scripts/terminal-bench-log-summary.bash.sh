IN=agent/inputs/logs/
OUT=agent/outputs/terminal-bench-log-summary.bash/
MANIFEST=${MANIFEST:-$PASH_TOP/evaluation/benchmarks/agent/logs.txt}
REF_DATE=${REF_DATE:-2025-08-12}

pure_func() {
    d="$1"
    grep -hoE '\[(ERROR|WARNING|INFO)\]' | sort | uniq -c |
        awk -v d="$d" '{ sev = $2; gsub(/[][]/, "", sev); print d "," sev "," $1 }'
}
export -f pure_func

for fname in $(cat "$MANIFEST"); do
    d="${fname%%_*}"                                   # date from YYYY-MM-DD_<src>.log
    python3 $PASH_TOP/aws/s3-get-object.py "$IN$fname" /dev/stdout | pure_func "$d" | python3 $PASH_TOP/aws/s3-put-object.py "${OUT}counts.${fname}.stdout" /dev/stdin
done

echo "Done"