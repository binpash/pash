#!/bin/bash

# 1.3: sort top first names
cat $IN | cut -d ' ' -f 1 | sort | uniq -c | sort -r
