IN=$PASH_TOP/evaluation/tests/input/1M.txt

mkfifo {1,2,3}grams

cat "$IN" |
  tr -cs A-Za-z '\n' |
  tr A-Z a-z |
  tee 3grams 2grams 1grams > /dev/null &

cat 1grams |
    sort |
    uniq -c |
    sort -rn > 1-grams.txt &

cat 2grams |
    sort |
    uniq -c |
    sort -rn > 2-grams.txt &

cat 3grams |
    sort |
    uniq -c |
    sort -rn # >> 3-grams.txt

rm {1,2,3}grams {1,2,3}-grams.txt
