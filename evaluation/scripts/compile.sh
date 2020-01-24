#!/bin/bash

# Find markdown files  in the current directory tree, compile  them to HTML, and
# serve them over the network

# Requires: pandoc

IN=./input/
OUT=./output/out.txt

find $IN -name '*.md' |    # Parallelizable, given a distributed FS
    xargs pandoc |         # xargs is higher-order, pandoc is third-party
    gzip > $OUT            # Compress the result
#   nc -l 80               # netcat could default-but-configurably parallelizable


