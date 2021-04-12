#!/bin/bash
# FIXME where is the zannotate binary zzzz
# tag: zannotate scan data
set -e
IN=${IN:-$PASH_TOP/evaluation/benchmarks/posh/input}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/posh/input/output}

export zannotate=/home/deeptir/go/bin/zannotate
#export filename=$IN/2019091303_port_80.json
export filename=$IN/scan_data/bash/test_80_40GB.json
export mrt_file=$IN/scan_data/2019-10-12.0500.mrt
export annotated=$IN/scan_data/bash/annotated
export as_popularity=$IN/scan_data/bash/as_popularity
export file1=$IN/scan_data/bash/file1
export file2=$IN/scan_data/bash/file2

cat $filename |  $zannotate -routing -routing-mrt-file=$mrt_file -input-file-type=json > $annotated
cat $annotated | jq ".ip" | tr -d '"' > $file1
cat $annotated | jq -c ".zannotate.routing.asn" > $file2
pr -mts, $file1 $file2 | awk -F',' "{ a[\$2]++; } END { for (n in a) print n \",\" a[n] } " | sort -k2 -n -t',' -r > $as_popularity
