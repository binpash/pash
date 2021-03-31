#!/bin/bash

if [[ "$1" == "-c" ]]; then
    rm -f genesis exodus
    exit
fi

# FIXME Could we apply these to all books across Gutenberg
if [ ! -f ./genesis ]; then
  curl -sf http://www.gutenberg.org/cache/epub/8001/pg8001.txt > genesis
fi 

if [ ! -f ./exodus ]; then
  curl -sf http://www.gutenberg.org/cache/epub/8001/pg8001.txt > exodus
fi
