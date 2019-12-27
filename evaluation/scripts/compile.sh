#!/bin/bash

# Find markdown files  in the current directory tree, compile  them to HTML, and
# serve them over the network

# Requires: pandoc

find . -name '*.md' |    # Parallelizable, given a distributed FS
    xargs pandoc |       # xargs is higher-order, trivially parallelizable; pandoc is third-party
    nc -l 80             # netcat could default-but-configurably parallelizable


