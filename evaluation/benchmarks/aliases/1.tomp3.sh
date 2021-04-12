#!/bin/bash

# tag: wav-to-mp3
set -e

IN=${WAV:-$PASH_TOP/evaluation/benchmarks/aliases/input/wav}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/aliases/input/out}

find $IN -name '*.wav' | 
    xargs -n1 basename |
    sed "s;\(.*\);-i $IN/\1 -ab 192000 $OUT/\1.mp3;" |
    xargs -L1  ffmpeg -y -loglevel quiet -hide_banner
