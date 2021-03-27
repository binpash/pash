#!/bin/bash

# Find all 2-grams in a piece of text
# FIXME: does not calculate frequencies

cat $IN |
  tr -cs A-Za-z '\n' |
  tr A-Z a-z |
  alt_bigrams_aux

