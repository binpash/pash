# Count vowel Sequences
INPUT=${INPUT:-$PASH_TOP/evaluation/scripts/input/genesis}
ls $IN/ | xargs cat | tr 'a-z' '[A-Z]' | tr -sc 'AEIOU' '[\012*]'| sort | uniq -c 
