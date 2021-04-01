# Count morphs in genesis
IN=${IN:-$PASH_TOP/evaluation/benchmarks/poets/input/pg/}
ls $IN | xargs cat | spell | sed 's/ .*//g' | sort | uniq -c 
