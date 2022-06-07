#!/bin/bash
#hadoop jar /homes/das160/hadoop-3.2.3/share/hadoop/tools/lib/hadoop-streaming-3.2.3.jar -files 6.sh -D mapred.reduce.tasks=0 -input /3.txt -output /unix506 -mapper 6.sh
cat $1 | cut -d ' ' -f 2 | cut -c 1-1 | tr -d '\n' | tr '[A-Z]' '[a-z]'

