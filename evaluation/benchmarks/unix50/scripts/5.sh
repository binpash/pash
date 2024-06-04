#!/bin/bash

# 2.1: get all Unix utilities
cat $IN | cut -d ' ' -f 4 | tr -d ','
