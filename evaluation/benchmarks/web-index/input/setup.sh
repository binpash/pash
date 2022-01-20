#!/bin/bash

set -e

cd $(dirname $0)

wiki_archive="https://dumps.wikimedia.org/other/static_html_dumps/current/en/wikipedia-en-html.tar.7z"
PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}

. "$PASH_TOP/scripts/utils.sh"

[ "$1" = "-c" ] && rm-files en/ *.7z *.tar  100.txt 1000.txt full small

if [ -e ./en/ ]; then
    exit 0
fi

if [ "$1" = "--small" ]; then
    # 100 entries
	wget http://pac-n4.csail.mit.edu:81/pash_data/small/web-index.small.zip
    unzip web-index.small.zip
    mv small/* .
    rm -rf small
elif [ "$1" = "--full" ]; then
    rm -rf en
    # 1000 entries
	wget http://pac-n4.csail.mit.edu:81/pash_data/full/web-index.full.zip
    unzip web-index.full.zip
    mv full/* .
    rm -rf full
elif [ "$1" = "--gen-full" ]; then
    wget $wiki_archive || eexit "cannot fetch wikipedia"
    7za x wikipedia-en-html.tar.7z
    tar -xvf wikipedia-en-html.tar
    wget http://ndr.md/data/wikipedia/index.txt || eexit "cannot fetch wiki indices"
fi
