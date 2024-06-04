#!/bin/bash

# 7.1: identify number of AT&T unix versions
cat $IN | cut -f 1 | grep 'AT&T' | wc -l
