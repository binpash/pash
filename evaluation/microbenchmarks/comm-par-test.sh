#!/bin/bash

# Tests the parallelization of comm with a configuration input and a stream input.
mkfifo s1 s2

## TODO: Add a grep after the second cat triggers an error
## TODO: Add a new line in the second s2 to have difference in output
## TODO: In theory we could parallelize a node that has a sequence in it!

cat $IN > s1 &
cat $IN > s2 &
comm -23 s1 s2

rm s1 s2
