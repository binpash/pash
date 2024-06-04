#!/bin/bash

# 1.2: extract names and sort
cat $IN | head -n 2 | cut -d ' ' -f 2
