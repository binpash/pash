#!/bin/bash
#hadoop jar /homes/das160/hadoop-3.2.3/share/hadoop/tools/lib/hadoop-streaming-3.2.3.jar -D mapred.reduce.tasks=0 -input /PaSh/9.4.txt -output /PaSh/unix527 -file 26.sh -mapper 26.sh
cat $1 | tr ' ' '\n' | grep "\"" | sed 4d | cut -d "\"" -f 2 | tr -d '\n'
