#!/bin/bash
#hadoop jar /homes/das160/hadoop-3.2.3/share/hadoop/tools/lib/hadoop-streaming-3.2.3.jar -files 7_map.sh,7_reduce.sh -input /4.txt -output /unix5072 -mapper 7_map.sh -reducer 7_reduce.sh
cat $1 | tr ' ' '\n' | grep '\.' | wc -l
