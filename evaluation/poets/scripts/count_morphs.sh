# Count morphs in genesis
INPUT=${INPUT:-$PATH_TOP/evaluation/scripts/input/genesis}
spell -v ${INPUT} | sed 's/ .*//g' | sort | uniq -c 
