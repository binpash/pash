#!/bin/bash

# 9.2: extract the word BELL
cat $IN | cut -c 1-1 | tr -d '\n'
