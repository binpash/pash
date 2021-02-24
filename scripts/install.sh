#!/usr/bin/env bash

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
#sudo add-apt-repository ppa:avsm/ppa
#sudo apt update
#sudo apt install -y opam m4

# Or just download opam
# wget # https://github.com/ocaml/opam/releases/download/2.0.3/opam-2.0.3-x86_64-linux
#opam init -y --compiler=4.07.1
#eval $(opam env)
#opam switch create 4.07.1
#eval $(opam env)
# ocaml -version | grep -o -E '[0-9]+.[0-9]+.[0-9]+$'

# Python 3.8 for older versions of Ubuntu
# sudo add-apt-repository ppa:deadsnakes/ppa

# FIXME At times, this will be 
# 
git submodule init
git submodule update

## If option -p is set, also run the sudo
if [ "$prepare_sudo_install_flag" -eq 1 ]; then
    echo "Running preparation sudo apt install:"
    echo "|-- running apt update..."
    sudo apt-get update &> $LOG_DIR/apt_update.log
    echo "|-- running apt install..."
    sudo apt-get install -y libtool m4 automake pkg-config libffi-dev python3 python3-pip wamerican-insane bc bsdmainutils &> $LOG_DIR/apt_install.log
else
    echo "Requires libtool, m4, automake, opam, pkg-config, libffi-dev, python3, pip for python3, a dictionary, bc, bsdmainutils"
    echo "Ensure that you have them by running:"
    echo "  sudo apt install libtool m4 automake pkg-config libffi-dev python3 python3-pip wamerican-insane bc bsdmainutils"
    echo -n "Press 'y' if you have these dependencies installed. "
    while : ; do
        read -n 1 k <&1
        if [[ $k = y ]] ; then
            echo ""
            echo "Proceeding..."
            break
        fi
    done
fi

## Old installation also required `sudo apt install opam` and `yes | opam init`
# sudo apt-get install opam
# yes | opam init &> $LOG_DIR/opam_init.log

# Build the parser (requires libtool, m4, automake, opam)
echo "Building parser..."
cd compiler/parser
echo "|-- making libdash... (requires sudo)"
## TODO: How can we get rid of that `sudo make install` in here?
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

# Install python3 dependencies
# 16.04 requires distutils, but has no python3-distutils
# sudo apt install python-distutils-extra
# sudo apt install python3-distutils-extra
# sudo apt install python3-distutils
# sudo apt remove python3-pip
# sudo python3 -m easy_install pip
echo "Installing python dependencies..."
python3 -m pip install jsonpickle &> $LOG_DIR/pip_install_jsonpickle.log
python3 -m pip install -U PyYAML &> $LOG_DIR/pip_install_pyyaml.log
python3 -m pip install numpy &> $LOG_DIR/pip_install_numpy.log
python3 -m pip install matplotlib &> $LOG_DIR/pip_install_matplotlib.log

# Generate inputs
echo "Generating input files..."
cd evaluation/scripts/input
./gen.sh
cd ../../../

# Export necessary environment variables
export PASH_TOP=$PWD

## This is necessary for the parser to link to libdash
echo "Do not forget to export LD_LIBRARY_PATH as shown below :)"
set -v
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/"
