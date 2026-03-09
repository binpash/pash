  sed "s#^#$HOME/wikipedia/#" |
  xargs cat |
  iconv -c -t ascii//TRANSLIT |
  pandoc +RTS -K64m -RTS --from html --to plain --quiet |
  tr -cs A-Za-z '\n' |
  tr A-Z a-z |
  grep -vwFf ../evaluation/scripts/web-search/stopwords.txt |
  ../evaluation/scripts/web-search/stem-words.js
