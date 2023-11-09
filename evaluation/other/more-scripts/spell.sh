#!/bin/bash

# Spell-check one or more man-pages
# (TODO: could use markdown pages?)
# 
# TODO: `groff is an interesting "pure", whose wrapper only needs split input
# TODO: files carefully.

# Data: atlas-group.cs.brown.edu/data/dummy/ronn.1
# dict depends on the system (and has to be sorted), so we assume it exists
dict=./input/dict.txt

# mercurial page
IN=/usr/share/man/man1/hg.1.gz
OUT=./output/out.txt

cat $IN |
  gunzip |                       # extract
  groff -t -e -mandoc -Tascii |  # remove formatting commands
  # prepare filename | 
  col -bx |                      # remove backspaces / linefeeds
  tr A-Z a-z |                   # map upper to lower case
  tr -d '[:punct:]' |            # remove punctuation
  sort |                         # put words in alphabetical order
  uniq |                         # remove duplicate words
  comm -13 $dict - > $OUT        # report words not in dictionary 

# The 2nd `tr` is also an example  of writing a stage in multiple ways
# some trivially parallelizable (tr), others not so much (awk)
# https://stackoverflow.com/questions/49199136/removing-punctuation-from-txt-file-in-linux-terminal-using-tr-and-awk-comman
