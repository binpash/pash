#FIXME: combine to a single regex?
# verses with 2 or more, 3 or more, exactly 2 instances of light.
INPUT=${INPUT:-$PASH_TOP/evaluation/scripts/input/genesis}
ls $IN/ | xargs cat | grep -c 'light.*light'  
ls $IN/ | xargs cat | grep -c 'light.*light.*light' 
ls $IN/ | xargs cat | grep 'light.*light' | grep -vc 'light.*light.*light'
#grep -c 'light.*light' ${INPUT} | paste -sd+ | bc
#grep -c 'light.*light.*light' ${INPUT} | paste -sd+ | bc
#grep 'light.*light' ${INPUT} | grep -vc 'light.*light.*light'  | paste -sd+ | bc
