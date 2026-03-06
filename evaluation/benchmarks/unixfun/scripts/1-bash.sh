#!/bin/bash

IN=unixfun/inputs/1_20G.txt
OUT=unixfun/outputs/
python3 aws/s3-get-object.py "$IN" /dev/stdout | cut -d ' ' -f 2 | sort | python3 aws/s3-put-object.py ${OUT}stdout.txt /dev/stdin