#!/bin/bash

PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}

[[ "$1" == "-c" ]] && { rm -rf genesis exodus pg; exit; }

setup_dataset() {
  if [ ! -f ./book_links.txt ]; then
    wget -q -O book_links.txt "https://atlas-group.cs.brown.edu/data/gutenberg/books.txt"
    if [ ! -f book_links.txt ]; then
        echo "Failed to download book_links.txt"
        exit 1
    fi
  fi

  if [ ! -f ./genesis ]; then
      curl -sf https://atlas-group.cs.brown.edu/data/gutenberg/8/0/0/8001/8001.txt > genesis
      "$PASH_TOP/scripts/append_nl_if_not.sh" genesis
  fi 

  if [ ! -f ./exodus ]; then
    curl -sf https://atlas-group.cs.brown.edu/data/gutenberg/3/3/4/2/33420/33420-0.txt > exodus
    "$PASH_TOP/scripts/append_nl_if_not.sh" exodus
  fi

  if [ ! -e ./pg ]; then
    mkdir pg
    cd pg
  if [[ "$1" == "--full" ]]; then
    cat ../book_links.txt | while IFS= read -r line
    do
        full_url="https://atlas-group.cs.brown.edu/data/gutenberg/${line}"
        echo "Downloading $full_url"
        wget -q "$full_url"
    done
  else
    book_count=40
    head -n $book_count ../book_links.txt | while IFS= read -r line
    do
        full_url="https://atlas-group.cs.brown.edu/data/gutenberg/${line}"
        echo "Downloading $full_url"
        wget -q "$full_url"
    done
  fi
    for f in *.txt; do
      "$PASH_TOP/scripts/append_nl_if_not.sh" $f
    done
    cd ..
  fi
}

source_var() {
  if [[ "$1" == "--small" ]]; then
    export ENTRIES=40
  else
    # 1% of the input
    export ENTRIES=1060
  fi
}
