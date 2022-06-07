#!/bin/bash
cat $1 | tr ' ' '\n' | grep 1969 | wc -l
