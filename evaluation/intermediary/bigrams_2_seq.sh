#!/bin/bash

# Find all 2-grams in a piece of text
# FIXME: does not calculate frequencies

cat $IN $IN |
  tr -cs A-Za-z '\n' |
  tr A-Z a-z |
  bigrams_aux |
  sort |
  uniq


