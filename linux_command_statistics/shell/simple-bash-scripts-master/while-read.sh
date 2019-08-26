#!/bin/bash
# while-read: read lines from a file
count=0
while read; do
	printf "%d %s\n" $REPLY
	count=$(expr $count + 1)
done <$1
