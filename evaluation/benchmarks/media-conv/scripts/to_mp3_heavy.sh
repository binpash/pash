#!/bin/bash
# tag: wav-to-mp3
# IN=${IN:-/media-conv/wav}
# OUT=${OUT:-$DISH_TOP/evaluation/media-conv/outputs/mp3}

pure_func(){
    ffmpeg -y -i pipe:0 -f mp3 -ab 192000 pipe:1  2>/dev/null
}
export -f pure_func

for item in $(seq 1 $ENTRIES); do
    for j in $(cat "$PASH_TOP/evaluation/benchmarks/media-conv/heavy_wav_list");do
        cat ${IN}${j} | pure_func > ${OUT}${j}.${item}.stdout.log; 
    done
done
echo 'done';