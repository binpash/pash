#!/bin/bash
IN=${IN:-/max-temp/temperatures.txt}

## Processing
hdfs dfs -cat "${IN}" |
  cut -c 89-92 |
  grep -v 999 |
  sort -rn |
  head -n1 > max.txt

hdfs dfs -cat "${IN}" |
  cut -c 89-92 |
  grep -v 999 |
  sort -n |
  head -n1 > min.txt

hdfs dfs -cat "${IN}" |
  cut -c 89-92 |
  grep -v 999 |
  awk "{ total += \$1; count++ } END { print total/count }" > average.txt 
