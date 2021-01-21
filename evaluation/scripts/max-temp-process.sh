#!/bin/bash

## Processing
cat $IN |
    cut -c 89-92 |
    grep -v 999 |
    sort -rn |
    head -n1
