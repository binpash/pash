#!/bin/bash
PROXY=$([ "$(hostname)" == "deathstar" ] && echo "gamma.ndr.md" || echo "localhost")
WIKI="$HOME/wikipedia/"
export WIKI
# Squash all HTML for each URL into a single line, streaming fashion
# It also prefixes with the URL

page_per_line () {
  cat "$WIKI/$0" | tr -d "\n\r" | tr -d '\n' | sed -e '/.$/a\'
}

export -f page_per_line

# xargs:
# add `-t` for debugging
cat $WIKI/index_h_100.txt | xargs -0 -d '\n' -n 1 bash -c 'page_per_line "$@"'
