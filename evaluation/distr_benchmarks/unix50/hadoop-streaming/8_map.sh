#!/bin/bash
#hadoop jar /homes/das160/hadoop-3.2.3/share/hadoop/tools/lib/hadoop-streaming-3.2.3.jar -files 8_map.sh,8_reduce.sh -input /4.txt -output /unix508 -mapper 8_map.sh -reducer 8_reduce.sh
cat $1 | tr ' ' '\n' | grep 'x' | grep '\.' | wc -l
