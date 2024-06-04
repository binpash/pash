#!/bin/bash

# 9.4: four corners with E centered, for an "X" configuration
cat $IN | tr ' ' '\n' | grep "\"" | sed 4d | cut -d "\"" -f 2 | tr -d '\n'
