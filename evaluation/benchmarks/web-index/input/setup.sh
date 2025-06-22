#!/bin/bash

#set -e

wiki_archive="https://dumps.wikimedia.org/other/static_html_dumps/current/en/wikipedia-en-html.tar.7z"
PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}

. "$PASH_TOP/scripts/utils.sh"

[ "$1" = "-c" ] && rm-files en/ *.7z *.tar  500.txt 1000.txt full small

setup_dataset() {
  rm -rf ../1-grams.txt ../2-grams.txt 
  
  ## Downloading the dataset needs to happen for both small and large
  if [[ ! -d ./en ]]; then
    wget $wiki_archive || eexit "cannot fetch wikipedia"
    7za x wikipedia-en-html.tar.7z
    tar -xvf wikipedia-en-html.tar
    # wget http://ndr.md/data/wikipedia/index.txt # || eexit "cannot fetch wiki indices"
    # It is actually OK if we don't have this index since we download the 500/1000 below
  fi

  if [ "$1" = "--small" ]; then
    # 500 entries
    # wget http://pac-n4.csail.mit.edu:81/pash_data/small/web-index.small.zip
    # unzip web-index.small.zip
    # mv small/500.txt .
    # rm -rf small web-index.small.zip

    wget https://atlas.cs.brown.edu/data/wikipedia/input_small/index.txt -O 500.txt || eexit "cannot fetch small index"
  else
    # elif [ "$1" = "--full" ]; then
    # the default full 
    # 1000 entries
    # wget http://pac-n4.csail.mit.edu:81/pash_data/full/web-index.full.zip
    # unzip web-index.full.zip
    # mv full/1000.txt .
    # rm -rf full web-index.full.zip
    wget "https://atlas.cs.brown.edu/data/wikipedia/index1g.txt" -O 1000.txt || eexit "cannot fetch full index"
  fi
}

source_var() {
  export WEB_INDEX_DIR=$PASH_TOP/evaluation/benchmarks/web-index/input
  export WIKI=$PASH_TOP/evaluation/benchmarks/web-index/input/
  if [[ "$1" == "--small" ]]; then
    export IN=$PASH_TOP/evaluation/benchmarks/web-index/input/500.txt
  else
    export IN=$PASH_TOP/evaluation/benchmarks/web-index/input/1000.txt
  fi
}
