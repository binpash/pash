#!/bin/bash
cat $1 | grep 'print' | cut -d "\"" -f 2 | cut -c 1-12
