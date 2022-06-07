#!/bin/bash
#hadoop jar /homes/das160/hadoop-3.2.3/share/hadoop/tools/lib/hadoop-streaming-3.2.3.jar -files 10_map.sh,10_reduce.sh -input /4.txt -output /unix510 -mapper 10_map.sh -reducer 10_reduce.sh
cat $1 | tr ' ' '\n' | grep '\.' | cut -d '.' -f 2 | cut -c 1-1 | tr '[a-z]' 'P'
