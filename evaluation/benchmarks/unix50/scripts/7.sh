#!/bin/bash

# 4.1: find number of rounds
cat $IN | tr ' ' '\n' | grep '\.' | wc -l
