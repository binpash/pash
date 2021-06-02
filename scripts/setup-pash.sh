#!/usr/bin/env bash

set -e

cd $(dirname $0)
PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}
cd $PASH_TOP

LOG_DIR=$PWD/install_logs
REQ=$PWD/scripts/requirements.txt
mkdir -p $LOG_DIR
git submodule init
git submodule update

echo "Building parser..."
cd compiler/parser
echo "|-- making libdash..."
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
virtualenv -p python3 pashenv
source pashenv/bin/activate
# after we activate the environment
pip install -r ${REQ} &> $LOG_DIR/pip_reqs.log

echo "Generating input files..."
$PASH_TOP/evaluation/tests/input/setup.sh

# export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/"
echo " * * * "
echo "Do not forget to export PASH_TOP before using pash: \`export PASH_TOP=$PASH_TOP\`"
echo '(optionally, you can update PATH to include it: `export PATH=$PATH:$PASH_TOP`)'
echo "Do not forget to run \`source pashenv/bin/activate\`!"
echo "Do not forget to run \`deactivate\` when leaving the pash environment!"
