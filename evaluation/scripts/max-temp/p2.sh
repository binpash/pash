#!/bin/bash

grep gz |
    # head -n 1 | # <-- remove this live for full scale
    tr -s ' ' |
    cut -d ' ' -f9

