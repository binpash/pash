#!/bin/bash

# 9.7: Four corners
cat $IN | sed 2d | sed 2d | tr -c '[A-Z]' '\n' | tr -d '\n'
