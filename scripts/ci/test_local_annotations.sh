#!/bin/bash

# Define the file path
FILE="$PASH_TOP/../annotations/pash_annotations/annotation_generation/AnnotationGeneration.py"

# Check if the file exists
if [[ ! -f "$FILE" ]]; then
    echo "Error: File not found at $FILE"
    exit 1
fi

# Comment out the line containing '"cat": "Cat"' if not already commented
sed -i 's/^\(\s*"cat": "Cat",\)/# \1/' "$FILE"

echo "Successfully commented out the 'cat' entry in $FILE"

