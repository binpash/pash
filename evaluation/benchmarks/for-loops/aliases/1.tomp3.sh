#!/bin/bash

# tag: wav-to-mp3
set -e
IN=${WAV:-$PASH_TOP/evaluation/benchmarks/aliases/input/wav}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/aliases/input/out}
# TODO: could we use something like this
rm -f $OUT/*
for FILE in $IN/*.wav ; 
do 
logname=$OUT/$(basename $FILE .wav).log
ffmpeg -y -i $FILE -f mp3 -ab 192000 $OUT/$(basename $FILE .wav).mp3 > '$OUT'/'$(basename $FILE .wav).log'; 
done
#find $IN -name '*.wav' | 
#    xargs -n1 basename |
#    sed "s;\(.*\);-i $IN/\1 -ab 192000 $OUT/\1.mp3;" |
#    xargs -L1  ffmpeg -y -loglevel quiet -hide_banner