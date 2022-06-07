#!/bin/bash
#hadoop jar /homes/das160/hadoop-3.2.3/share/hadoop/tools/lib/hadoop-streaming-3.2.3.jar -D mapred.reduce.tasks=0 -input /PaSh/9.9.txt -output /PaSh/unix531 -file 31.sh -mapper 31.sh
cat $1 | tr -c '[a-z][A-Z]' '\n' | grep '[A-Z]' | sed 1d | sed 1d | sed 2d | sed 3d | sed 5d | tr -c '[A-Z]' '\n' | tr -d '\n'
