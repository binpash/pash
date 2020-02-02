
cat "$IN_DIR/p1.out_2_00" "$IN_DIR/p1.out_2_01" |
  sed "s#^#$WIKI#" |
  xargs cat |
  iconv -c -t ascii//TRANSLIT |
  pandoc +RTS -K64m -RTS --from html --to plain --quiet |
  tr -cs A-Za-z '\n' |
  tr A-Z a-z |
  grep -vwFf $WEB_INDEX_DIR/stopwords.txt |
  $WEB_INDEX_DIR/stem-words.js 
