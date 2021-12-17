#!/bin/bash

set -e

cd $(dirname $0)

wiki_archive="https://dumps.wikimedia.org/other/static_html_dumps/current/en/wikipedia-en-html.tar.7z"
PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}

. "$PASH_TOP/scripts/utils.sh"

[ "$1" = "-c" ] && rm-files en/ *.7z *.tar

if [ ! -e ./en/ ]; then
  wget $wiki_archive || eexit "cannot fetch wikipedia"
  7za x wikipedia-en-html.tar.7z
  tar -xvf wikipedia-en-html.tar
  wget http://ndr.md/data/wikipedia/index.txt || eexit "cannot fetch wiki indices"
fi

