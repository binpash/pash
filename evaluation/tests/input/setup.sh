#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

# Source utils for helper functions
. "$REPO_ROOT/scripts/utils.sh"

# Always set PASH_TOP to the correct location (src/pash)
# This overrides any stale values from old PaSh installations
export PASH_TOP="$REPO_ROOT/src/pash"
cd "$SCRIPT_DIR"

distro=$(printf '%s\n' "$distro" | LC_ALL=C tr '[:upper:]' '[:lower:]')
head_sz="M"
# now do different things depending on distro
case "$distro" in
   freebsd*) 
    head_sz="m"
    ;;
esac

[ "$1" = "-c" ] && rm-files 1M.txt all_cmds.txt words sorted_words 10M.txt

if [ ! -f ./1M.txt ]; then
  curl -sf --connect-timeout 10 'atlas.cs.brown.edu/data/dummy/1M.txt' > 1M.txt
  if [ $? -ne 0 ]; then
    curl -f 'https://zenodo.org/record/7650885/files/1M.txt' > 1M.txt
    if [ $? -ne 0 ]; then
      curl -sf 'http://www.gutenberg.org/files/2600/2600-0.txt' | head -c 1${head_sz} > 1M.txt
      [ $? -ne 0 ] && eexit 'cannot find 1M.txt'
    fi
  fi
  append_nl_if_not ./1M.txt
fi

if [ ! -f ./all_cmds.txt ]; then
  if [ "$(hostname)" = "deathstar" ]; then
    curl -sf --connect-timeout 10 'atlas.cs.brown.edu/data/dummy/all_cmds.txt' > all_cmds.txt
    if [ $? -ne 0 ]; then
      curl -f 'https://zenodo.org/record/7650885/files/all_cmds.txt' > all_cmds.txt || eexit "all_cmds not found"
    fi
  else
    ls /usr/bin/* > all_cmds.txt
  fi
  append_nl_if_not ./all_cmds.txt
fi

if [ ! -f ./words ]; then
  curl -sf --connect-timeout 10 'atlas.cs.brown.edu/data/dummy/words' > words
  if [ $? -ne 0 ]; then
    curl -f 'https://zenodo.org/record/7650885/files/words' > words
    if [ $? -ne 0 ]; then
      if [ $(uname) = 'Darwin' ]; then
        cp /usr/share/dict/web2 words || eexit "cannot find dict file"
      else
        # apt install wamerican-insane
        cp /usr/share/dict/words words || eexit "cannot find dict file"
      fi
    fi
  fi
  append_nl_if_not words
fi

if [ ! -f ./10M.txt ]; then
  touch 10M.txt
  for (( i = 0; i < 10; i++ )); do
    cat 1M.txt >> 10M.txt
  done
fi

if [ ! -f ./ab.txt ]; then
  touch ab.txt
  for (( i = 0; i < 10; i++ )); do
    echo 'aaa' >> ab.txt
    echo 'bbb' >> ab.txt
  done
fi


## Re-sort words for this machine
if [ ! -f ./sorted_words ]; then
  sort words > sorted_words
fi
