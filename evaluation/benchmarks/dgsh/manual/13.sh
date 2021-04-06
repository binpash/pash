#!/bin/bash
# tag: venue author compare
# from: https://www2.dmst.aueb.gr/dds/sw/dgsh/#author-compare
set -e

IN=${IN:-$PASH_TOP/evaluation/benchmarks/dgsh/input/}
OUT=${OUT:-$PASH_TOP/evaluation/benchmarks/dgsh/input}
cd ${OUT}

# Extract and sort author names
sorted_authors()
{
  sed -n 's/<author>\([^<]*\)<\/author>/\1/p' |
  sort
}

# Escape a string to make it a valid sed(1) pattern
escape()
{
  echo "$1" | sed 's/\([/\\]\)/\\\1/g'
}

export -f sorted_authors

rm -f sed_a sed_b a b a_ch1 a_ch2 b_ch1 b_ch2 
mkfifo sed_a sed_b a b a_ch1 a_ch2 b_ch1 b_ch2 

#
# Output the two venue authors as two output streams
cat sed_a | sed   "
/^<.*key=\"$(escape $1)/,/<title>/ w >|" > a &

cat sed_b | sed  "
/^<.*key=\"$(escape $2)/,/<title>/ w >|" > b  &
cat ${IN}mini.xml| tail -n 1000 | tee sed_a sed_b  > /dev/null 
# at this point we have the results of seds, time to make more tee

# fork the sed outputs as inputs to the other 4 commands
cat a | tee a_ch1  > /dev/null &
cat b | tee b_ch1  > /dev/null &



echo -n "$1 papers: " &&  cat a_ch1 | grep -c '^<.* mdate=.* key=' | cat > a_ch2 &
echo -n "$2 papers: " &&  cat b_ch1 | grep -c '^<.* mdate=.* key=' | cat > b_ch2 &

# redirect sed to the sorted_authors
cat a_ch2 | sorted_authors  &
cat b_ch2 | sorted_authors  &
# Incompete
