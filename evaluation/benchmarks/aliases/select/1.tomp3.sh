#!/bin/bash

# tag: wav-to-mp3
set -e

IN=$PASH_TOP/evaluation/benchmarks/aliases/meta/wav
OUT=$PASH_TOP/evaluation/benchmarks/aliases/meta/out
find $IN -name '*.wav' | xargs -I {} -P 16 ffmpeg -y -hide_banner -loglevel error -i {} -ab 192000 $OUTPUT/{}.mp3
