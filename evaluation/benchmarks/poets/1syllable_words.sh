# 1-syllable words 
# FIXES: name exercise number
INPUT=${IN:-$PASH_TOP/evaluation/scripts/input/genesis}
ls $IN/ | xargs cat | tr -sc '[A-Z][a-z]' '[ 12*]' | grep -i '^[^aeiou]*[aeiou][^aeiou]*$' | sort | uniq -c | sed 5q
