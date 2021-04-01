#!/bin/bash

# tag: wav-to-mp3
set -e

IN=$PASH_TOP/evaluation/benchmarks/aliases/meta/wav
OUT=$PASH_TOP/evaluation/benchmarks/aliases/meta/out

#find $IN -name '*.wav' | 
#    xargs -n1 basename |
#    sed "s;\(.*\);-i $IN/\1 -ab 192000 $OUT/\1.mp3;" |
#    xargs -0 -n1 echo ffmpeg

find . -iname "*.wav" -printf "%p -ab 192000 ${OUT}/%f.mp3\n" | 
    xargs -r -n4 ffmpeg -i
