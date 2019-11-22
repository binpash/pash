#!/bin/bash

set -x

IN=../scripts/input/i1G.txt
rm -f "#file8"
rm -f "#file11"
rm -f "#file10"
rm -f "#file12"
rm -f "#file9"
mkfifo "#file8"
mkfifo "#file11"
mkfifo "#file10"
mkfifo "#file14"
mkfifo "#file13"
mkfifo "#file12"
mkfifo "#file9"
cat $IN > "#file9" &
cat $IN > "#file10" &
cat "#file9" | tr A-Z a-z > "#file11" &
cat "#file10" | tr A-Z a-z > "#file12" &
cat "#file11" | sort > "#file13" &
cat "#file12" | sort > "#file14" &
# wait
sort -m "#file13" "#file14" > "#file8" &
cat "#file8" > "/tmp/distr_output"
