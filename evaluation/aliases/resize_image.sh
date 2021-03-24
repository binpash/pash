# Resize image 
INPUT=${INPUT:-$PASH_TOP/evaluation/aliases/input/jpg}
OUTPUT=${OUTPUT:-$PASH_TOP/evaluation/aliases/output}
cd $INPUT
find . -name "*.jpg" -o -name "*.png" -o -name "*.jpeg" -o -name "*.JPG" -o -name "*.PNG" \
    -o -name "*.JPEG" | xargs -P 16 -I {} sh -c "convert -resize 70% {} {}.70; mv {}.70 $OUTPUT" sh {}
