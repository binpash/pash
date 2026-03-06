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
    for file in "$directory"/*.hash
    do
        # Copy the file to the hash folder
        cp "$file" "$hash_folder"
    done
fi

# Loop through all directories in the parent directory
for folder in "outputs"/*/
do
    # Remove trailing slash
    folder=${folder%/}

    echo "Verifying folder: $folder"

    # Loop through all .hash files in the current directory
    for file in "$folder"/*.hash
    do
        # Extract the filename without the directory path and extension
        filename=$(basename $file)

        # Compare the hash with the hash in the hashes directory
        if ! diff "$hash_folder/$filename" "$folder/$filename";
        then
            # Print the filename and hash if they don't match
            echo "File: $folder/$filename hash diff failed!"
        fi
    done
done
