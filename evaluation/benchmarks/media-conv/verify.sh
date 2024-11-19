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


    # Loop through all .out files in the directory and subdirectories
    find "$directory" -type f -name "*.hash" | while read -r file;
    do
        echo $file
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

    # Loop through all .hash files in the current directory and subdirectories
    find "$folder" -type f -name "*.hash" | while read -r file;
    do
        # Extract the filename without the directory path and extension
        filename=$(basename $file)

        # Compare the hash with the hash in the hashes directory
        if ! diff "$hash_folder/$filename" "$file";
        then
            # Print the filename and hash if they don't match
            echo "File: $file hash diff failed comparing to $hash_folder/$filename"
        fi
    done
done
