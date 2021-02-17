# Count morphs in genesis
INPUT=${INPUT:-$PASH_TOP/evaluation/scripts/input/poets/genesis}
spell -v ${INPUT} | sed 's/ .*//g' | sort | uniq -c 
