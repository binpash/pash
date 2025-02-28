# #!/bin/bash

# # Define the requirements file path
# REQ_FILE="requirements.txt"

# # Determine PaSh top directory
# PASH_TOP=$(dirname "$(realpath "$0")")

# # Default annotations directory (sibling of PASH_TOP)
# DEFAULT_ANNOTATIONS_DIR="$(dirname "$PASH_TOP")/annotations"

# # Use user-provided path if given, otherwise use default
# if [ -n "$2" ]; then
#     ANNOTATIONS_DIR="$2"
#     ANNOTATIONS_DIR=$(realpath "$ANNOTATIONS_DIR")
# else
#     ANNOTATIONS_DIR="$DEFAULT_ANNOTATIONS_DIR"
# fi

# echo "📂 Annotations directory set to: $ANNOTATIONS_DIR"

# # Check if at least one argument is provided
# if [ $# -lt 1 ]; then
#     echo "Usage: $0 <local|pypi> [annotations_dir]"
#     exit 1
# fi

# # Store the current version of pash-annotations only if using pypi mode or if it is for the first time
# if [ "$1" == "pypi" ] && [ ! -s pash-annotations.version ]; then
#     grep 'pash-annotations' requirements.txt > pash-annotations.version
# fi


# # Ensure the annotations directory exists for local mode
# if [ "$1" == "local" ]; then
#     if [ ! -d "$ANNOTATIONS_DIR" ]; then
#         echo "🔍 $ANNOTATIONS_DIR not found. Cloning annotations repository..."
#         git clone https://github.com/binpash/annotations.git "$ANNOTATIONS_DIR"

#         if [ $? -ne 0 ]; then
#             echo "❌ Failed to clone annotations repository."
#             exit 1
#         fi
#     else
#         echo "✅ $ANNOTATIONS_DIR already exists. Skipping clone."
#     fi

#     # Get the absolute path
#     ANNOTATIONS_PATH=$(realpath "$ANNOTATIONS_DIR")
# fi

# # Generate `requirements.txt` dynamically
# echo "🔄 Generating $REQ_FILE..."

# # Copy base requirements
# cp requirements.base.txt $REQ_FILE

# # If local mode is selected, add the editable annotation repo
# if [ "$1" == "local" ]; then
#     echo "-e $ANNOTATIONS_PATH" >> $REQ_FILE
#     echo "✅ Local annotations added to $REQ_FILE"
# else
#     # Ensure `pash-annotations` is only added once
#     if ! grep -q 'pash-annotations' $REQ_FILE; then
#         cat pash-annotations.version >> $REQ_FILE
#     fi
#     echo "✅ Using PyPI annotations: $(cat pash-annotations.version)"
# fi

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

