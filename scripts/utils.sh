#!/usr/bin/env bash

# #Check that we are in the appropriate directory where setup.sh is
# #https://stackoverflow.com/a/246128
# DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
# echo "changing to $DIR to run setup.sh"
# cd $DIR

# another solution for capturing HTTP status code
# https://superuser.com/a/590170

eexit(){
  echo $1 'please email pash-devs@googlegroups.com'
  exit 1
}

nargs(){
  echo $# $1 $2
}

rm-files(){
  echo "${@}"
  rm -r "${@}"
  exit 0
}

append_nl_if_not(){
  ## Adds a newline at the end of a file if it doesn't already end in a newline.
  ## Used to prepare inputs for PaSh.
  if [ -z "$1" ]; then
    echo "No file argument given!"
    exit 1
  else
    if [ ! -f "$1" ]; then
      echo "File $1 doesn't exist!"
      exit 1
    else
      tail -c 1 "$1" | od -ta | grep -q nl
      if [ $? -eq 1 ]; then
        echo >> "$1"
      fi
    fi
  fi
}

install_deps_source_setup() {
    # move to the input directory
    cd input/
    # check if there are dependencies
    if [ -e install-deps.sh ]; then
        echo "Installing dependencies"
        bash install-deps.sh
    fi
    # source the setup file
    # it contains the fetch dataset function
    # and the export variable function for IN, IN_PRE
    source setup.sh
    # fetch the dataset
    setup_dataset $1 > /dev/null
    cd ..
}
#########################
# The command line help #
#########################
usage() {
    echo "Usage: `basename $0` [option...] -- shell script to build PaSh"
    echo
    echo "   -h, --help                 Show this help message"
    echo "   -o, --opt-agg              Install g++-10 and switch to it as main compiler. Build the optimized c++ aggregators (run with sudo)"
    echo "   -s, --show-deps            Show all the required dependencies (does not setup/deploy PaSh nor its dependencies)"
    echo "   -e, --install-eval         Install all the dependencies needed for reproducing the evaluation figures (uses sudo, only for Ubuntu/Debian currently)"
    echo
    exit 1
}

##########################################
# Install all the required libraries and #
# dependencies for PaSh evaluation       #
##########################################
install_eval_deps() {
    echo "Installing evaluation dependencies (needs sudo)"
    # needed for majority of the benchmarks (not available in docker instances)
    sudo apt-get install unzip
    paths="$(find $PASH_TOP/evaluation/benchmarks -name install-deps.sh)"
    for f in $(echo $paths); do
        path=$(dirname $(readlink -f $f))
        cd $path
        bash install-deps.sh
        cd - > /dev/null
    done
    echo "Generating PDF plots of the evaluation results is optional and requires R-packages"
    echo "Follow Installation Guide from: $PASH_TOP/evaluation/eval_script/README.md"
}

##########################################
# parse and read the command line args   #
##########################################
read_cmd_args() {
    # Transform long options to short ones
    for arg in "$@"; do
        shift
        case "$arg" in
            "--opt-agg")        set -- "$@" "-o" ;;
            "--show-deps")      set -- "$@" "-s" ;;
            "--install-eval")   set -- "$@" "-e" ;;
            "--help")           set -- "$@" "-h" ;;
            *)                  set -- "$@" "$arg"
        esac
    done

    while getopts 'opsreh' opt; do
        case $opt in
            # passthrough the variable to the Makefile for libdash
            o) export optimized_agg_flag=1 ;;
            s) export show_deps=1 ;;
            r) export show_eval_deps=1 ;;
            e) export install_eval=1 ;;
            h) usage >&2 ;;
            *) echo 'Error in command line parsing' >&2
                exit 1
        esac
    done
}
