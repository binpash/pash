# generate thumbnails with size 100x100 
INPUT=${INPUT:-$PASH_TOP/evaluation/scripts/input/posh/jpg/}
OUTPUT=${OUTPUT:-$PASH_TOP/evaluation/scripts/input/posh/output}
find ${INPUT} -name "*.jpg" | xargs -n 1 mogrify -format jpg -path ${OUTPUT} -thumbnail 100x100
