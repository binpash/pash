#!/bin/bash
# tag: resize image 
# IN=${JPG:-/media-conv/jpg}
# OUT=${OUT:-$DISH_TOP/evaluation/media-conv/outputs/jpg}
# mkdir -p $2

pure_func () {
     convert -resize 70% "-" "-"
}
export -f pure_func

for input in $(cat "$PASH_TOP/evaluation/benchmarks/media-conv/jpg-list.txt" | head -n ${ENTRIES})
do
    cat $IN$input | pure_func > $OUT$input.out
done

echo 'done';
