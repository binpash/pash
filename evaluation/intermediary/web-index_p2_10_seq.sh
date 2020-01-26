
cat "$IN_DIR/p1.out_10_00" "$IN_DIR/p1.out_10_01" "$IN_DIR/p1.out_10_02" "$IN_DIR/p1.out_10_03" "$IN_DIR/p1.out_10_04" "$IN_DIR/p1.out_10_05" "$IN_DIR/p1.out_10_06" "$IN_DIR/p1.out_10_07" "$IN_DIR/p1.out_10_08" "$IN_DIR/p1.out_10_09" |
  sed "s#^#$WIKI#" |
  xargs cat |
  iconv -c -t ascii//TRANSLIT |
  pandoc --from html --to plain --quiet |
  tr -cs A-Za-z '\n' |
  tr A-Z a-z |
  grep -vwFf $WEB_INDEX_DIR/stopwords.txt |
  $WEB_INDEX_DIR/stem-words.js 
