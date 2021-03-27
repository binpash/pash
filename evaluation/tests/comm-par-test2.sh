#!/bin/bash

# Tests the parallelization of comm with a configuration input and a stream input.
mkfifo s1 s2

cat $IN > s1 &
cat $IN | grep "king" > s2 &
comm -23 - s2 < s1

rm s1 s2
