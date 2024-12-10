#!/usr/bin/env bash

for file in logs_original/*; do
    file=$(basename $file)
    echo "Comparing: $file"
    sdiff -l <(grep "caruca" logs_original/$file) <(grep "caruca" logs_caruca/$file)
done
