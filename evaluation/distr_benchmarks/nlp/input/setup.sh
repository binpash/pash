#!/bin/bash

PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}

[[ "$1" == "-c" ]] && { rm -rf genesis exodus pg; exit; }

if [ ! -f ./genesis ]; then
    curl -sf https://www.gutenberg.org/cache/epub/8001/pg8001.txt > genesis
    "$PASH_TOP/scripts/append_nl_if_not.sh" genesis
fi 

if [ ! -f ./exodus ]; then
  curl -sf https://www.gutenberg.org/files/33420/33420-0.txt > exodus
  "$PASH_TOP/scripts/append_nl_if_not.sh" exodus
fi

if [ ! -e ./pg ]; then
  mkdir pg
  cd pg
  if [[ "$1" == "--full" ]]; then
    echo 'N.b.: download/extraction will take about 10min'
    wget atlas-group.cs.brown.edu/data/pg.tar.xz # FIXME: moving to PG soon
    if [ $? -ne 0 ]; then
		cat <<-'EOF' | sed 's/^ *//'
		Downloading input dataset failed, thus need to manually rsync all books from  project gutenberg:
		rsync -av --del --prune-empty-dirs --include='*.txt' --include='*/' --exclude='*' ftp@ftp.ibiblio.org::gutenberg .
		please contact the pash developers pash-devs@googlegroups.com
		EOF
    exit 1
  fi
  cat pg.tar.xz | tar -xJ
  
  else
    # wget http://pac-n4.csail.mit.edu:81/pash_data/nlp.zip
    # unzip nlp.zip
    # mv data/* .
    # rm nlp.zip data -rf
    
    # Mock 1
    for (( i = 0; i < 60; i++ )); do
        touch "$i".txt
        cat ../genesis >> "$i".txt
    done
    # Mock 2
    for (( i = 61; i < 120; i++ )); do
        touch "$i".txt
        cat ../exodus >> "$i".txt
    done
  fi

  for f in *.txt; do
    "$PASH_TOP/scripts/append_nl_if_not.sh" $f
  done
  cd ..
  
fi

# Put files in hdfs
hdfs dfs -mkdir /nlp
hdfs dfs -put exodus /nlp/exodus
hdfs dfs -put genesis /nlp/genesis
hdfs dfs -put pg /nlp
