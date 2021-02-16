# Recursive PaSh calls to trigrams (defined below)
# trigram file scripts ???
INPUT=${INPUT:-$PATH_TOP/evaluation/scripts/input/genesis}
grep 'the land of' ${INPUT} | sh trigram | sort -nr | sed 5q 
grep 'And he said' ${INPUT} | sh trigram | sort -nr | sed 5q
