#!/bin/bash

# Ensure that the script fails if something failed
set -e

LOG_DIR=$PWD/install_logs
mkdir -p $LOG_DIR


prepare_sudo_install_flag=0
while getopts 'p' opt; do
    case $opt in
        p) prepare_sudo_install_flag=1 ;;
        *) echo 'Error in command line parsing' >&2
           exit 1
    esac
done
shift "$(( OPTIND - 1 ))"

## for earlier versions of Debian/Ubuntu
## https://github.com/janestreet/install-ocaml
#add-apt-repository ppa:avsm/ppa
#apt update
#apt install -y opam m4
#opam init -y --compiler=4.07.1
#eval $(opam env)
# ocaml -version | grep -o -E '[0-9]+.[0-9]+.[0-9]+$'

git submodule init
git submodule update

## If option -p is set, also run the sudo
if [ "$prepare_sudo_install_flag" -eq 1 ]; then
    echo "Running preparation sudo apt install and opam init:"
    echo "|-- running apt update..."
    sudo apt-get update &> $LOG_DIR/apt_update.log
    echo "|-- running apt install..."
    sudo apt-get install -y libtool m4 automake opam pkg-config libffi-dev python3.8 python3-pip &> $LOG_DIR/apt_install.log
    opam -y init &> $LOG_DIR/opam_init.log
    # opam update
else
    echo "Requires libtool, m4, automake, opam, pkg-config, libffi-dev, python3.8, pip for python3"
    echo "Ensure that you have them by running:"
    echo "  sudo apt install libtool m4 automake opam pkg-config libffi-dev python3.8 python3-pip"
    echo "  opam init"
    echo "Press 'y' if you have these dependencies installed."
    while : ; do
        read -n 1 k <&1
        if [[ $k = y ]] ; then
            echo "Proceeding..."
            break
        fi
    done
fi


# Build the parser (requires libtool, m4, automake, opam)
echo "Building parser..."
eval $(opam config env)
cd comipler/parser
echo "|-- installing opam dependencies..."
make opam-dependencies &> $LOG_DIR/make_opam_dependencies.log
echo "|-- making libdash..."
make libdash &> $LOG_DIR/make_libdash.log
echo "|-- making parser..."
make &> $LOG_DIR/make.log
cd ../../

echo "Building runtime..."
# Build runtime tools: eager, split
cd runtime/
make &> $LOG_DIR/make.log
cd ../

# Install python3.8 dependencies
echo "Installing python dependencies..."
python3.8 -m pip install jsonpickle &> $LOG_DIR/pip_install_jsonpickle.log
python3.8 -m pip install -U PyYAML &> $LOG_DIR/pip_install_pyyaml.log

# Generate inputs
echo "Generating input files..."
cd evaluation/scripts/input
./gen.sh
cd ../../../

# Export necessary environment variables
export PASH_TOP=$PWD
export PASH_PARSER=${PASH_TOP}/parser/parse_to_json.native

## This is necessary for the parser to link to libdash
echo "Do not forget to export LD_LIBRARY_PATH as shown below :)"
set -v
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/"
