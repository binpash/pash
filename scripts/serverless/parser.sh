# for file in tmp folder
# read the file, grep the 169.254.76.1:XXXX line

#!/bin/bash

# # Directory containing files
# DIR="tmp-32"

# # Iterate over each file in the tmp folder
# for file in "$DIR"/*; do
#     # Check if it is a regular file
#     if [[ -f "$file" ]]; then
#         echo "Processing file: $file"
        
#         # Read the file and grep for lines with 169.254.76.1:XXXX pattern
#         grep -E "169\.254\.76\.1:[0-9]{4,5}" "$file"
#     fi
# done

#!/bin/bash

# Directory containing files
DIR="tmp"

# Iterate over each file in the tmp folder
for file in "$DIR"/*; do
    # Check if it is a regular file
    if [[ -f "$file" ]]; then

        
        # Grep for the exact pattern and only output the matched portion
        grep -o "169\.254\.76\.1:[0-9]\{4,5\}" "$file"
    fi
done