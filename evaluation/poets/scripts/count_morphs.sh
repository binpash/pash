# Count morphs in genesis
INPUT=${INPUT:-inputs/genesis}
spell -v ${INPUT} | sed 's/ .*//g' | sort | uniq -c 
