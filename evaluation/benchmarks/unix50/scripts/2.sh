#!/bin/bash

# 1.1: extract names and sort
cat $IN | cut -d ' ' -f 2 | sort >${OUT}stdout.txt