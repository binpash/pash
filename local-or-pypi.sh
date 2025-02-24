#!/bin/bash

# Define the requirements file path
REQ_FILE="requirements.txt"
ANNOTATIONS_DIR="$HOME/annotations"

# Check if an argument is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <local|pypi>"
    exit 1
fi

# If using local mode, ensure ~/annotations exists
if [ "$1" == "local" ]; then
    if [ ! -d "$ANNOTATIONS_DIR" ]; then
        echo "🔍 ~/annotations not found. Cloning annotations repository..."
        git clone https://github.com/binpash/annotations.git "$ANNOTATIONS_DIR"

        # Check if cloning was successful
        if [ $? -ne 0 ]; then
            echo "❌ Failed to clone annotations repository."
            exit 1
        fi
    else
        echo "✅ ~/annotations already exists. Skipping clone."
    fi

    # Get the absolute path
    ANNOTATIONS_PATH=$(realpath "$ANNOTATIONS_DIR")
fi

# Write the correct `requirements.txt`
if [ "$1" == "pypi" ]; then
    echo "📦 Writing PyPI requirements to $REQ_FILE..."
    cat <<EOF > "$REQ_FILE"
graphviz
libdash
pash-annotations==0.2.2
shasta==0.1.0
sh-expand>=0.1.6
EOF
elif [ "$1" == "local" ]; then
    echo "💻 Writing local development requirements to $REQ_FILE..."
    cat <<EOF > "$REQ_FILE"
graphviz
libdash
shasta==0.1.0
sh-expand>=0.1.6
-e $ANNOTATIONS_PATH
EOF
else
    echo "❌ Invalid argument: $1"
    echo "Usage: $0 <local|pypi>"
    exit 1
fi

echo "✅ Requirements file updated successfully."

