#!/bin/bash
#hadoop jar <some prefix>hadoop-streaming-<vers>.jar -files spell_map.sh,spell_reduce.sh -input <path to input file in HDFS> -output <output-dir in HDFS> -mapper spell_map.sh -reducer spell_reduce.sh
cat $1 | iconv -f utf-8 -t ascii//translit | col -bx | tr -c 'A-Za-z' '[\n*]' | grep -v "^\s*$"  | tr A-Z a-z | tr -d '[:punct:]'
