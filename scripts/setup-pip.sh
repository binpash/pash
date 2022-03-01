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
# save distro in the init file
echo "export distro=$distro" > ~/.pash_init
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
esac

cd ../

# echo "Generating input files..."
# $PASH_TOP/evaluation/tests/input/setup.sh
cp /opt/pash/pa.sh /usr/bin/

echo " * * * "
echo "Do not forget to export PASH_TOP before using pash: \`export PASH_TOP=$PASH_TOP\`"
echo '(optionally, you can update PATH to include it: `export PATH=$PATH:$PASH_TOP`)'
echo " * * * "