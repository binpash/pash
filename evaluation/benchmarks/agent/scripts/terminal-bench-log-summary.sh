IN=${IN:-$PASH_TOP/evaluation/benchmarks/agent/inputs/logs/}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/agent/outputs/}
MANIFEST=${MANIFEST:-$PASH_TOP/evaluation/benchmarks/agent/logs.txt}
REF_DATE=${REF_DATE:-2025-08-12}

mkdir -p "$OUT"

pure_func() {
    d="$1"
    grep -hoE '\[(ERROR|WARNING|INFO)\]' | sort | uniq -c |
        awk -v d="$d" '{ sev = $2; gsub(/[][]/, "", sev); print d "," sev "," $1 }'
}
export -f pure_func

for fname in $(cat "$MANIFEST"); do
    d="${fname%%_*}"                                   # date from YYYY-MM-DD_<src>.log
    cat "${IN}${fname}" | pure_func "$d" > "${OUT}counts.${fname}.stdout"
done

echo "Done"