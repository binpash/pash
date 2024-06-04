#!/bin/bash

# 5.1: extract hello world
cat $IN | grep 'print' | cut -d "\"" -f 2 | cut -c 1-12
