#!/bin/bash

#set -e

PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}

. "$PASH_TOP/scripts/utils.sh"

[ "$1" = "-c" ] && rm-files en/ *.7z *.tar  500.txt 1000.txt full small

setup_dataset() {
  rm -rf ../1-grams.txt ../2-grams.txt

  ## Downloading the dataset needs to happen for both small and large
  if [[ -d ./articles ]]; then
    exit 0
  fi

  if [ "$1" = "--small" ]; then
    wiki_archive="https://atlas.cs.brown.edu/data/wikipedia/wikipedia1g.tar.gz"
    wget -O wikipedia.tar.gz $wiki_archive || eexit "cannot fetch wikipedia"
    tar -xf wikipedia.tar.gz -C .
    mv articles1g articles
    wget -O index.txt https://atlas.cs.brown.edu/data/wikipedia/index1g.txt
  fi

  if [ "$1" = "--full" ]; then
    wiki_archive="https://atlas.cs.brown.edu/data/wikipedia/wikipedia10g.tar.gz"
    wget $wiki_archive -O wikipedia10g.tar.gz || eexit "cannot fetch full wiki"
    tar -xvzf wikipedia10g.tar.gz
    rm wikipedia10g.tar.gz
    wget -O index.txt https://atlas.cs.brown.edu/data/wikipedia/index10g.txt
  fi

}

source_var() {
  export WEB_INDEX_DIR=$PASH_TOP/evaluation/benchmarks/web-index/input
  export WIKI=$PASH_TOP/evaluation/benchmarks/web-index/input/articles/
  if [[ "$1" == "--small" ]]; then
    export IN=$PASH_TOP/evaluation/benchmarks/web-index/input/index.txt
  else
    export IN=$PASH_TOP/evaluation/benchmarks/web-index/input/index.txt
  fi
}