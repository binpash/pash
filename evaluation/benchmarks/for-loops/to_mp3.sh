#!/bin/bash
# tag: wav-to-mp3
set -e
IN=${IN:-$PASH_TOP/evaluation/benchmarks/for-loops/input/wav}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/for-loops/input/output/mp3}
LOGS=${OUT}/logs
mkdir -p ${LOGS}

for FILE in $IN/*.wav ; 
do 
    ffmpeg -y -i $FILE -f mp3 -ab 192000 $OUT/$(basename $FILE).mp3 > $LOGS/$(basename $FILE).log; 
done
