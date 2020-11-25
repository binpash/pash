#!/bin/bash

# Ensure that the script fails if something failed
set -e

echo "Requires libtool, m4, automake, opam, pkg-config, libffi-dev, python3.8"
echo "Ensure that you have them by running:"
echo "  sudo apt install libtool m4 automake opam pkg-config libffi-dev python3.8"
echo "  opam init"
echo "Press 'y' if you have these dependencies installed."
while : ; do
    read -n 1 k <&1
    if [[ $k = y ]] ; then
        echo "Proceeding..."
        break
    fi
done

# Build the parser (requires libtool, m4, automake, opam)
echo "Building parser..."
eval $(opam config env)
cd parser
make opam-dependencies
make libdash
make
cd ../

echo "Building runtime..."
# Build runtime tools: eager, split
cd evaluation/tools/
make
cd ../../

# Install python3.8 dependencies
python3.8 -m pip install jsonpickle
python3.8 -m pip install -U PyYAML

# Generate inputs
echo "Generating input files"
cd evaluation/scripts/input
./gen.sh
cd ../../../

# Export necessary environment variables
export PASH_TOP=$PWD
export PASH_PARSER=${PASH_TOP}/parser/parse_to_json.native