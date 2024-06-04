#!/bin/bash

# 3.1: get lowercase first letter of last names (awk)
cat $IN | cut -d ' ' -f 2 | cut -c 1-1 | tr -d '\n' | tr '[A-Z]' '[a-z]'
