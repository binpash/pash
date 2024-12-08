#!/usr/bin/env bash

for file in logs_original/*; do
    file=$(basename $file)
    echo "Comparing: $file"
    diff <(grep "caruca class" logs_original/$file) <(grep "caruca class" logs_caruca/$file)
done
