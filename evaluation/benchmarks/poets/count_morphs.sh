# Count morphs in genesis
INPUT=${INPUT:-$PASH_TOP/evaluation/scripts/input/genesis}
ls $IN/ | xargs cat | spell | sed 's/ .*//g' | sort | uniq -c 
