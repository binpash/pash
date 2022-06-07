#!/bin/bash
jarpath="/opt/hadoop-3.2.2/share/hadoop/tools/lib/hadoop-streaming-3.2.2.jar"
basepath=""

# Bigrams
# Diff
hadoop jar $jarpath -files nfa-regex.sh -D mapred.reduce.tasks=0 -input $basepath/1G.txt -output /nfa-regex -mapper nfa-regex.sh
# Set-diff
hadoop jar $jarpath -files shortest-scripts_map.sh,shortest-scripts_reduce.sh -input $basepath/1G.txt -output /shortest-scripts -mapper shortest-scripts_map.sh -reducer shortest-scripts_reduce.sh
hadoop jar $jarpath -files sort-sort_map.sh,sort-sort_reduce.sh -input $basepath/1G.txt -output /sort-sort -mapper sort-sort_map.sh -reducer sort-sort_reduce.sh
hadoop jar $jarpath -files sort.sh -D mapred.reduce.tasks=0 -input $basepath/1G.txt -output /sort -mapper sort.sh
hadoop jar $jarpath -files spell_map.sh,spell_reduce.sh -input $basepath/1G.txt -output /spell -mapper spell_map.sh -reducer spell_reduce.sh
hadoop jar $jarpath -files top-n_map.sh,top-n_reduce.sh -input $basepath/1G.txt -output /top-n -mapper top-n_map.sh -reducer top-n_reduce.sh
hadoop jar $jarpath -files wf_map.sh,wf_reduce.sh -input $basepath/1G.txt -output /wf -mapper wf_map.sh -reducer wf_reduce.sh
