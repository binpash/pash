#!/bin/bash
#!/bin/bash
# source: posh benchmarks

BENCHMARK_DIR="$(cd "$(dirname "$0")/.." && pwd)"

filename=${1:-$IN}
[ -d "$filename" ] && filename="$filename/all_logs.jsonl"

mrt_file=${2:-${MRT:-"$BENCHMARK_DIR/inputs/routeviews.mrt"}}

if [ $# -ge 3 ]; then
    annotated="$3"; file1="$4"; file2="$5"; as_popularity="$6"
else
    mkdir -p "$OUT"
    annotated="$OUT/annotated.jsonl"
    file1="$OUT/file1"
    file2="$OUT/file2"
    as_popularity="$OUT/as_popularity.csv"
fi

cat "$filename" | zannotate -routing -routing-mrt-file="$mrt_file" -input-file-type=json \
  | jq -r '"\(.ip),\(.zannotate.routing.asn)"' \
  | cut -d',' -f2 \
  | sort \
  | uniq -c \
  | awk '{print $2","$1}' \
  | sort -k2 -n -t',' -r \
  > "$as_popularity"

# cat "$filename" | zannotate -routing -routing-mrt-file=$mrt_file -input-file-type=json > "$annotated"
# cat "$annotated" | jq ".ip" | tr -d '"' > "$file1"
# cat "$annotated" | jq -c ".zannotate.routing.asn" > "$file2"
# pr -mts, $file1 $file2 | awk -F',' "{ a[\$2]++; } END { for (n in a) print n \",\" a[n] } " | sort -k2 -n -t',' -r > "$as_popularity"
