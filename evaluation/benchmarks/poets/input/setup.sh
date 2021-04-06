#!/bin/bash

PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}

mkdir -p $PASH_TOP/evaluation/benchmarks/poets/input/output
[[ "$1" == "-c" ]] && { rm -rf genesis exodus pg; exit; }

if [ ! -f ./genesis ]; then
  curl -sf http://www.gutenberg.org/cache/epub/8001/pg8001.txt > genesis
  "$PASH_TOP/scripts/append_nl_if_not.sh" genesis
fi 

if [ ! -f ./exodus ]; then
  curl -sf http://www.gutenberg.org/cache/epub/8001/pg8001.txt > exodus
  "$PASH_TOP/scripts/append_nl_if_not.sh" exodus
fi

if [ ! -e ./pg ]; then
  mkdir pg
  cd pg
  echo 'N.b.: download/extraction will take about 10min'
  curl -f ndr.md/data/pg.tar.xz > pg.tar.xz
  if [ $? -ne 0 ]; then
    cat <<-'EOF'
    Downloading input dataset failed, thus need to manually rsync all books from  project gutenberg:
    rsync -av --del --prune-empty-dirs --include='*.txt' --include='*/' --exclude='*' ftp@ftp.ibiblio.org::gutenberg .
    please contact the pash developers pash-devs@googlegroups.com
    EOF
  fi
  cat pg.tar.gz | tar -xJ
  for f in *.txt; do
    "$PASH_TOP/scripts/append_nl_if_not.sh" $f
  done
  cd ..
fi

