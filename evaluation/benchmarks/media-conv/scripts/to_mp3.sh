#!/bin/bash
# tag: wav-to-mp3
# IN=${IN:-/media-conv/wav}
# OUT=${OUT:-$DISH_TOP/evaluation/media-conv/outputs/mp3}

pure_func(){
    ffmpeg -y -i pipe:0 -f mp3 -ab 192000 pipe:1  2>/dev/null
}
export -f pure_func

for item in $(cat "$PASH_TOP/evaluation/benchmarks/media-conv/wav-list.txt" | head -n ${ENTRIES});
do
    cat $IN$item | pure_func > $OUT$item.out
done

echo 'done';
