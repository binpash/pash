
## TODO: For now we have to output everything to stdout and not to
## different files when distributing (1-grams.txt)

cat '1grams' |
    sort |
    uniq -c |
    sort -rn &

cat '2grams' |
    tr -cs A-Za-z '\n' |
    tr A-Z a-z |
    wi_bigrams_aux |
    uniq -c |
    sort -rn &

cat '3grams' |
    tr -cs A-Za-z '\n' |
    tr A-Z a-z |
    wi_trigrams_aux |
    uniq -c |
    sort -rn &

# s2---NLP manipulation:  get text
cat "$IN_DIR/web-index_p1.out" |
  pandoc --from html --to plain --quiet |
  tr -cs A-Za-z '\n' |
  tr A-Z a-z |
  iconv -c -t ascii//TRANSLIT |
  grep -vwFf stopwords.txt |
  ./stem-words.js | # stem-to-roots
  tee '3grams' '2grams' '1grams' > /dev/null &

wait
