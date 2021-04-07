#!/bin/bash
# tag: word_properties
# from: https://www2.dmst.aueb.gr/dds/sw/dgsh/#word-properties
set -e
# FIXME Output is weird
IN=${MINI:-$PASH_TOP/evaluation/benchmarks/dgsh/input/mini.xml}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/dgsh/input}
cd ${OUT}
#OUTPUT=${OUTPUT:-$PASH_TOP/evaluation/benchmarks/dgsh/output}
rm -f a b c d a1 b1 c1 d1 
mkfifo a b c d a1 b1 c1 d1 
cat a | cat > a1 &
# List two-letter palindromes
cat b | sed 's/.*\(.\)\(.\)\2\1.*/p: \1\2-\2\1/;t
    g' | cat  > b1 &

# List four consecutive consonants
cat c | sed -E 's/.*([^aeiouyAEIOUY]{4}).*/c: \1/;t
    g' | cat > c1 &

# List length of words longer than 12 characters
cat d | awk '{if (length($1) > 12) print "l:", length($1);
    else print ""}' | cat > d1 &

cat ${IN} | tr -cs a-zA-Z \\n | sort -u | tee a b c d  > /dev/null  
paste a1 b1 c1 d1  | fgrep :
#cat d1 | fgrep :
rm a b c d a1 b1 c1 d1
