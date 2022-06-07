#!/bin/bash

# hadoop jar /homes/das160/hadoop-3.2.3/share/hadoop/tools/lib/hadoop-streaming-3.2.3.jar -input /1.txt -output /unix502 -file 1.sh -mapper 1.sh
cat $1 | cut -d ' ' -f 2
