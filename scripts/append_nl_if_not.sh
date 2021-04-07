#!/bin/bash

## Adds a newline at the end of a file if it doesn't already end in a newline.
## Used to prepare inputs for PaSh.

if [ -z "$1" ]; then
    echo "No file argument given!"
    exit 1
else
    if [ ! -f "$1" ]; then
        echo "File $1 doesn't exist!"
        exit 1
    else
        tail -c 1 "$1" | od -ta | grep -q nl
        if [ $? -eq 1 ]
        then
        echo >> "$1"
        fi
    fi
fi
