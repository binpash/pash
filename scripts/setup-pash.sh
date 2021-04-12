#!/usr/bin/env bash

set -e

cd $(dirname $0)
PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}
cd $PASH_TOP

LOG_DIR=$PWD/install_logs
mkdir -p $LOG_DIR

git submodule init
git submodule update

echo "Building parser..."
cd compiler/parser
echo "|-- making libdash... (requires sudo)"
make libdash &> $LOG_DIR/make_libdash.log
cd ../../

## This was the old parser installation that required opam.
# # Build the parser (requires libtool, m4, automake, opam)
# echo "Building parser..."
# eval $(opam config env)
# cd compiler/parser
# echo "|-- installing opam dependencies..."
# make opam-dependencies &> $LOG_DIR/make_opam_dependencies.log
# echo "|-- making libdash... (requires sudo)"
# ## TODO: How can we get rid of that `sudo make install` in here?
# make libdash &> $LOG_DIR/make_libdash.log
# make libdash-ocaml &>> $LOG_DIR/make_libdash.log
# echo "|-- making parser..."
# make &> $LOG_DIR/make.log
# cd ../../


echo "Building runtime..."
# Build runtime tools: eager, split
cd runtime/
make &> $LOG_DIR/make.log
cd ../

echo "Installing python dependencies..."
python3 -m pip install jsonpickle &> $LOG_DIR/pip_install_jsonpickle.log
python3 -m pip install -U PyYAML &> $LOG_DIR/pip_install_pyyaml.log
python3 -m pip install numpy &> $LOG_DIR/pip_install_numpy.log
python3 -m pip install matplotlib &> $LOG_DIR/pip_install_matplotlib.log

echo "Generating input files..."
$PASH_TOP/evaluation/tests/input/setup.sh

## This is necessary for the parser to link to libdash
echo "Do not forget to export PASH_TOP as shown below :)"
set -v
# export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/"
export PASH_TOP=$PASH_TOP

