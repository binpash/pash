#!/bin/bash

cd "$(realpath $(dirname "$0"))"
mkdir -p inputs
cd inputs

inputs=(1 3 10 11 12 2 4 5 6 7 8 9.1 9.2 9.3 9.4 9.5 9.6 9.7 9.8 9.9 10 11)

# inputs=(9.9)

S3_BUCKET_PREFIX="s3://$AWS_BUCKET"
S3_BENCHMARK_DIR="unix50"
S3_INPUTS_DIR="$S3_BENCHMARK_DIR/inputs"

for input in ${inputs[@]}
do

    echo "Processing ${input}.txt"

    if [ ! -f "${input}.txt" ]; then
        wget "http://atlas-group.cs.brown.edu/data/unix50/${input}.txt" -q
    fi

    if [ ! -f "${input}_1M.txt" ]; then
        file_content_size=$(wc -c < "${input}.txt")
        iteration_limit=$((1048576 / $file_content_size))
        echo "Iteration limit: $iteration_limit, File content size: $file_content_size"
        for (( i = 0; i < iteration_limit; i++ )); do
            cat "${input}.txt" >> "${input}_1M.txt"
        done
    fi

    if [ ! -f "${input}_20G.txt" ]; then
        for (( i = 0; i < 20000; i++ )); do
            cat "${input}_1M.txt" >> "${input}_20G.txt"
        done
    fi

    echo "Uploading ${input}_20G.txt"
    aws s3 cp "${input}_20G.txt" "${S3_BUCKET_PREFIX}/${S3_INPUTS_DIR}/${input}_20G.txt"

    echo "Remove ${input}_20G.txt"
    rm "${input}_20G.txt"

    echo "Finished processing ${input}.txt"
done
