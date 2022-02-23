#!/usr/bin/env bash

set -e

cd $(dirname $0)
# check the git status of the project
if git rev-parse --git-dir > /dev/null 2>&1; then
    # we have cloned from the git repo, so all the .git related files/metadata are available
    git submodule init
    git submodule update
    # set PASH_TOP
    PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}
else
    # set PASH_TOP to the root folder of the project if it is not available
    PASH_TOP=${PASH_TOP:-$PWD/..}
    # remove previous installation if it exists
    rm -rf $PASH_TOP/compiler/parser/libdash
    # we are in package mode, no .git information is available
    git clone https://github.com/angelhof/libdash/ $PASH_TOP/compiler/parser/libdash
fi
cd $PASH_TOP

LOG_DIR=$PASH_TOP/install_logs
mkdir -p $LOG_DIR
PYTHON_PKG_DIR=$PASH_TOP/python_pkgs
# remove the folder in case it exists
rm -rf $PYTHON_PKG_DIR
# create the new folder
mkdir -p $PYTHON_PKG_DIR

echo "Building parser..."
cd compiler/parser

if type lsb_release >/dev/null 2>&1 ; then
    distro=$(lsb_release -i -s)
elif [ -e /etc/os-release ] ; then
    distro=$(awk -F= '$1 == "ID" {print $2}' /etc/os-release)
fi

echo "|-- making libdash..."
# convert to lowercase
distro=$(printf '%s\n' "$distro" | LC_ALL=C tr '[:upper:]' '[:lower:]')
# now do different things depending on distro
case "$distro" in
    freebsd*) 
        gsed -i 's/ make/ gmake/g' Makefile
        gmake libdash &> $LOG_DIR/make_libdash.log
        echo "Building runtime..."
        # Build runtime tools: eager, split
        cd ../../runtime/
        gmake &> $LOG_DIR/make.log
        ;;
    *)
        make libdash &> $LOG_DIR/make_libdash.log
        echo "Building runtime..."
        # Build runtime tools: eager, split
        cd ../../runtime/
        make &> $LOG_DIR/make.log
        ;;
esac

# save distro in the init file
echo "export distro=$distro" > ~/.pash_init

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

cd ../

echo "Installing python dependencies..."
python3 -m pip install jsonpickle --root $PYTHON_PKG_DIR --force-reinstall #&> $LOG_DIR/pip_install_jsonpickle.log
python3 -m pip install pexpect --root $PYTHON_PKG_DIR --force-reinstall #&> $LOG_DIR/pip_install_pexpect.log
python3 -m pip install graphviz --root $PYTHON_PKG_DIR --force-reinstall #&> $LOG_DIR/pip_install_jsonpickle.log
python3 -m pip install numpy #&> $LOG_DIR/pip_install_numpy.log
python3 -m pip install matplotlib #&> $LOG_DIR/pip_install_matplotlib.log
# clean the python packages
cd $PYTHON_PKG_DIR
find . -name "site-packages" -type d | xargs -I '{}' mv $PYTHON_PKG_DIR/{}/ .
mv site-packages/* . && rm -rf site-packages

echo "Generating input files..."
$PASH_TOP/evaluation/tests/input/setup.sh

# export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/"
echo " * * * "
echo "Do not forget to export PASH_TOP before using pash: \`export PASH_TOP=$PASH_TOP\`"
echo '(optionally, you can update PATH to include it: `export PATH=$PATH:$PASH_TOP`)'

