#!/bin/bash

# 10.1: count Turing award recipients while working at Bell Labs
cat $IN | sed 1d | grep 'Bell' | cut -f 2 | wc -l
