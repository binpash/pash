#!/bin/bash
IN=../scripts/input/i1M.txt
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
cat "#file11" | sort  > "#file13" &
cat "#file12" | sort  > "#file14" &
sort -m "#file13" "#file14" > "#file8" &
cat "#file8" > "/tmp/distr_output"
