#!/bin/bash
# tag: compression_bench.sh
# from: https://www.spinellis.gr/sw/dgsh/#compress-compare
set -e
# IN wget https://dblp.uni-trier.de/xml/dblp.xml.gz
IN=${IN:-$PASH_TOP/evaluation/benchmarks/dgsh/input/}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/dgsh/input/}
cd ${OUT}
rm -f a b c d
mkfifo a b c d
cat a | wc -c | sed 's/^/Original size:\t/' &
cat b | xz -c | wc -c | sed 's/^/xz:\t\t/' &
cat c | bzip2 -c | wc -c | sed 's/^/bzip2:\t\t/' &
cat d | gzip -c | wc -c | sed 's/^/gzip:\t\t/' &
cat ${IN}mini.xml |tail -n 100 | tee a b c d > /dev/null
rm a b c d
