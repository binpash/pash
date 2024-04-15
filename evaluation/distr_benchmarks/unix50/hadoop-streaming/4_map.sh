#!/bin/bash

# hadoop jar /opt/hadoop-3.4.0/share/hadoop/tools/lib/hadoop-streaming-3.4.0.jar -files 4_map.sh,4_reduce.sh -input /1.txt -output /unix504 -mapper 4_map.sh -reducer 4_reduce.sh
cat $1 | cut -d ' ' -f 1
