#!/bin/bash

# Usage: ./setup-annotations.sh [ANNOTATIONS_DIR] [-f|--force]

FORCE=0
USER_PROVIDED_DIR=""

for arg in "$@"; do
  if [ "$arg" = "-f" ] || [ "$arg" = "--force" ]; then
    FORCE=1
  else
    USER_PROVIDED_DIR="$arg"
  fi
done

# Determine annotations directory
if [ -n "$USER_PROVIDED_DIR" ]; then
    ANNOTATIONS_DIR=$(realpath "$USER_PROVIDED_DIR")
else
    if [ -z "$PASH_TOP" ]; then
        echo "Error: PASH_TOP is not set and no annotations directory was provided."
        echo "Set PASH_TOP or provide a directory path as an argument."
        exit 1
    fi
    # PASH_TOP points to src/pash, so repo root is ../..
    REPO_ROOT="$(cd "$PASH_TOP/../.." && pwd)"
    ANNOTATIONS_DIR="$REPO_ROOT/annotations"
fi

echo "Annotations directory set to: $ANNOTATIONS_DIR"

# Check if directory exists
if [ -d "$ANNOTATIONS_DIR" ]; then
    if [ "$FORCE" -eq 1 ]; then
        echo "Removing existing annotations directory due to --force..."
        rm -rf "$ANNOTATIONS_DIR"
    else
        echo "Error: Annotations directory already exists at $ANNOTATIONS_DIR"
        echo "Use --force or -f to overwrite it."
        exit 1
    fi
fi

# Clone the annotations repository
echo "Downloading annotations repository..."
git clone https://github.com/binpash/annotations.git "$ANNOTATIONS_DIR"

if [ $? -ne 0 ]; then
    echo "Failed to clone annotations repository."
    exit 1
fi

ANNOTATIONS_PATH=$(realpath "$ANNOTATIONS_DIR")
echo "Annotations repository cloned to: $ANNOTATIONS_PATH"
