#!/bin/bash

pure_func() {
    tr -d '\n' | cut -c 1-4
}

# 9.1: extract the word PORT
cat $IN | tr ' ' '\n' | grep '[A-Z]' | tr '[a-z]' '\n' | grep '[A-Z]' | pure_func >${OUT}stdout.txt
