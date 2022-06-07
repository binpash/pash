#!/bin/bash
# hadoop jar /opt/hadoop-3.2.2/share/hadoop/tools/lib/hadoop-streaming-3.2.2.jar -files 5.sh -D mapred.reduce.tasks=0 -input /2.txt -output /unix506 -mapper 5.sh
cat $1 | cut -d ' ' -f 4 | tr -d ','
