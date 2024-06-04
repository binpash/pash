#!/bin/bash

# 9.6: Follow the directions for grep
cat $IN | tr ' ' '\n' | grep '[A-Z]' | sed 1d | sed 3d | sed 3d | tr '[a-z]' '\n' | grep '[A-Z]' | sed 3d | tr -c '[A-Z]' '\n' | tr -d '\n'
