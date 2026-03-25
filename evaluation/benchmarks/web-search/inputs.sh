#!/bin/bash

BENCHMARK_DIR=$(cd "$(dirname "$0")" && pwd)
input_dir=$BENCHMARK_DIR/inputs
URL=https://atlas.cs.brown.edu/data
mkdir -p "$input_dir"

suffix=""
is_small=false
is_min=false
for arg in "$@"; do
    case "$arg" in
        --small) is_small=true; suffix="_small" ;;
        --min)   is_min=true;   suffix="_min"   ;;
    esac
done

if [ ! -f "$input_dir/stopwords.txt" ]; then
    wget --no-check-certificate "${URL}/web-index/stopwords.txt" -O "$input_dir/stopwords.txt"
fi

if [ -d "$input_dir/articles$suffix" ]; then
    echo "Dataset already exists."
    exit 0
fi

if $is_min; then
    echo "Downloading the min dataset."
    wget --no-check-certificate -O "$input_dir/wikipedia$suffix.tar.gz" "${URL}/wikipedia/input_small/articles.tar.gz"
    wget --no-check-certificate -O "$input_dir/index$suffix.txt" "${URL}/wikipedia/input_small/index.txt"
    echo "Extracting the min dataset."
    tar -xf "$input_dir/wikipedia$suffix.tar.gz" -C "$input_dir" --no-same-owner
    mv "$input_dir/articles" "$input_dir/articles$suffix"
elif $is_small; then
    echo "Downloading the small dataset."
    wget --no-check-certificate -O "$input_dir/wikipedia$suffix.tar.gz" "${URL}/wikipedia/wikipedia1g.tar.gz"
    wget --no-check-certificate -O "$input_dir/index$suffix.txt" "${URL}/wikipedia/index1g.txt"
    echo "Extracting the small dataset."
    tar -xf "$input_dir/wikipedia$suffix.tar.gz" -C "$input_dir" --no-same-owner
    mv "$input_dir/articles1g" "$input_dir/articles$suffix"
else
    echo "Downloading the full dataset."
    wget --no-check-certificate -O "$input_dir/wikipedia$suffix.tar.gz" "${URL}/wikipedia/wikipedia10g.tar.gz"
    wget --no-check-certificate -O "$input_dir/index$suffix.txt" "${URL}/wikipedia/index10g.txt"
    echo "Extracting the full dataset."
    tar -xf "$input_dir/wikipedia$suffix.tar.gz" -C "$input_dir" --no-same-owner
    mv "$input_dir/articles10g" "$input_dir/articles$suffix"
fi
rm "$input_dir/wikipedia$suffix.tar.gz"
