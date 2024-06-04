#!/bin/bash

# 8.1: count unix birth-year
cat $IN | tr ' ' '\n' | grep 1969 | wc -l
