#!/usr/bin/env bash

# Ensure that the script fails if something failed
set -e

cd $(dirname $0)
PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}

LOG_DIR=$PASH_TOP/install_logs
mkdir -p $LOG_DIR

prepare_sudo_install_flag=0

#########################
# The command line help #
#########################
usage() {
    echo "Usage: `basename $0` [option...] -- shell script to build PaSh"
    echo
    echo "   -h, --help                 Show this help message"
    echo "   -o, --opt-agg              Install g++-10 and switch to it as main compiler. Build the optimized c++ aggregators (uses sudo)"
    echo "   -p, --prepare              Install all PaSh dependencies to the system (uses sudo)"
    echo "   -s, --show-deps            Show all the required dependencies (does not setup/deploy PaSh nor its dependencies)"
    echo "   -e, --install-eval         Install all the dependencies needed for reproducing the evaluation figures (uses sudo, only for Ubuntu/Debian currently)"
    echo
    # echo some stuff here for the -a or --add-options 
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

# Transform long options to short ones
for arg in "$@"; do
  shift
  case "$arg" in
    "--prepare")        set -- "$@" "-p" ;;
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
    p) prepare_sudo_install_flag=1 ;;
    s) export show_deps=1 ;;
    r) export show_eval_deps=1 ;;
    e) export install_eval=1 ;;
    h) usage >&2 ;;
    *) echo 'Error in command line parsing' >&2
      exit 1
  esac
done
shift "$(( OPTIND - 1 ))"

if [ "$show_deps" = 1 ]; then
    # just trigger the distro-deps script
    # no package is installed
    echo "Required packages for PaSh:" 
    ./distro-deps.sh
    exit 0
fi

if [ "$install_eval" -eq 1 ]; then
    install_eval_deps
fi
## If option -p is set, also run the sudo
if [ "$prepare_sudo_install_flag" -eq 1 ]; then
  # if we are within docker, we do not need to run sudo
  if [[ $PASH_HOST == "docker" ]]; then
    ./distro-deps.sh $optimized_agg_flag
  else
    # this is the default option for native installation
    sudo ./distro-deps.sh $optimized_agg_flag
  fi
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

./setup-pash.sh
