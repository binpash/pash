#!/bin/bash
#hadoop jar /homes/das160/hadoop-3.2.3/share/hadoop/tools/lib/hadoop-streaming-3.2.3.jar -files 9_map.sh,9_reduce.sh -input /4.txt -output /unix509 -mapper 9_map.sh -reducer 9_reduce.sh
cat $1 | tr ' ' '\n' | grep 'x' | grep '\.' | cut -d '.' -f 2 | grep -v '[KQRBN]' | wc -l
