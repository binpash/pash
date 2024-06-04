#!/bin/bash

# 9.1: extract the word PORT
cat $IN | tr ' ' '\n' | grep '[A-Z]' | tr '[a-z]' '\n' | grep '[A-Z]' | tr -d '\n' | cut -c 1-4
