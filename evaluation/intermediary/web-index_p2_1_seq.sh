
cat "$IN_DIR/p1.out" |
  sed "s#^#$WIKI#" |
  xargs cat |
  iconv -c -t ascii//TRANSLIT |
  pandoc --from html --to plain --quiet |
  tr -cs A-Za-z '\n' |
  tr A-Z a-z |
  grep -vwFf $WEB_INDEX_DIR/stopwords.txt |
  $WEB_INDEX_DIR/stem-words.js 
