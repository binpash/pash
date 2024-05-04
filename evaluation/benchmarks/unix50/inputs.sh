#!/bin/bash

cd "$(dirname "$0")" || exit 1

inputs=(
    1 2 3 4 5 6 7 8 9.1 9.2 9.3 9.4 9.5 9.6 9.7 9.8 9.9 10 11 12
)

input_dir="inputs"

mkdir -p $input_dir

for input in "${inputs[@]}"
do
    if [ ! -f "${input}.txt" ]; then
        wget "http://atlas-group.cs.brown.edu/data/unix50/${input}.txt" -q -P $input_dir
    fi

    cd $input_dir || exit 1

    # Concatenate file with itself until it reaches 1GB
    while [ "$(du -m "${input}.txt" | cut -f1)" -lt 1024 ]; do
        cat "${input}.txt" "${input}.txt" >> "${input}_temp.txt"
        mv -f "${input}_temp.txt" "${input}.txt"
    done

    # If the file size exceeds 1GB, split it into 1GB chunks
    if [ "$(du -m "${input}.txt" | cut -f1)" -gt 1025 ]; then
        split -C 1024m --numeric-suffixes "${input}.txt" "${input}_"
        mv -f "${input}_00" "${input}.txt"
        rm "${input}_"*
    fi

    cd ..

    echo "Finished processing ${input}.txt"
done
