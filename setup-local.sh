#!/bin/bash

# Define the requirements file path
REQ_FILE="requirements.txt"

# Determine PaSh top directory
PASH_TOP=$(dirname "$(realpath "$0")")

# Default annotations directory (sibling of PASH_TOP)
DEFAULT_ANNOTATIONS_DIR="$(dirname "$PASH_TOP")/annotations"

# Use user-provided path if given, otherwise use default
if [ -n "$1" ]; then
    ANNOTATIONS_DIR="$1"
    ANNOTATIONS_DIR=$(realpath "$ANNOTATIONS_DIR")
else
    ANNOTATIONS_DIR="$DEFAULT_ANNOTATIONS_DIR"
fi

echo "📂 Annotations directory set to: $ANNOTATIONS_DIR"

# Always clone the repository to the specified or default location
echo "🔄 Downloading annotations repository..."
rm -rf "$ANNOTATIONS_DIR"
git clone https://github.com/binpash/annotations.git "$ANNOTATIONS_DIR" 

if [ $? -ne 0 ]; then
    echo "❌ Failed to clone annotations repository."
    exit 1
fi

# Get the absolute path
ANNOTATIONS_PATH=$(realpath "$ANNOTATIONS_DIR")
echo $ANNOTATIONS_PATH > local-annotations-path.txt
echo "✅ Annotations repository cloned to: $ANNOTATIONS_PATH"
grep -qxF "$ANNOTATIONS_PATH" requirements.txt || echo -e "\n-e $ANNOTATIONS_PATH" >> requirements.txt

