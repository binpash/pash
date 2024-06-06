#!/bin/bash

# Exit immediately if a command exits with a non-zero status
# set -e

cd "$(realpath $(dirname "$0"))"

mkdir -p hashes/small

if [[ "$@" == *"--small"* ]]; then
    hash_folder="hashes/small"
else
    hash_folder="hashes"
fi

if [[ "$@" == *"--generate"* ]]; then
    # Directory to iterate over
    directory="outputs/bash"

    # Loop through all .out files in the directory
    find "$directory" -mindepth 2 -type f -name '*.out' | while read -r file;
    do
        # Extract the filename and dirname
        filename=$(basename "$file" .out)
        dirname=$(dirname "${file#$directory/}")

        # Generate SHA-256 hash
        hash=$(shasum -a 256 "$file" | awk '{ print $1 }')

        # Create subdirectory if not already
        mkdir -p $hash_folder/$dirname

        # Save the hash to a file
        echo "$hash" > "$hash_folder/$dirname/$filename.hash"

        # Print the filename and hash
        echo "File: $hash_folder/$dirname/$filename.hash | SHA-256 Hash: $hash"
    done
fi

# Loop through all directories in the parent directory
for folder in "outputs"/*
do
    echo "Verifying folder: $folder"

    # Loop through all .out files in the current directory
    find "$folder" -mindepth 2 -type f -name '*.out' | while read -r file;
    do
        # Extract the filename and dirname
        filename=$(basename "$file" .out)
        dirname=$(basename "$(dirname "$file")") # is the script_name

        if [ ! -f "$folder/$dirname/$filename.hash" ]; then
            # Generate SHA-256 hash
            hash=$(shasum -a 256 "$file" | awk '{ print $1 }')

            # Save the hash to a file
            echo "$hash" > "$folder/$dirname/$filename.hash"
        fi

        # Compare the hash with the hash in the hashes directory
        diff "$hash_folder/$dirname/$filename.hash" "$folder/$dirname/$filename.hash"

        # Print the filename and hash
        echo "File: $dirname/$filename | SHA-256 Hash: $(cat "$folder/$dirname/$filename.hash")"
    done
done
