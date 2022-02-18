#!/usr/bin/env bash

# Ensure that the script fails if something failed
set -e

cd $(dirname $0)
PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}

LOG_DIR=$PASH_TOP/install_logs
mkdir -p $LOG_DIR

prepare_sudo_install_flag=0

# Transform long options to short ones
for arg in "$@"; do
  shift
  case "$arg" in
    "--prepare")    set -- "$@" "-p" ;;
    "--opt-agg")    set -- "$@" "-o" ;;
    *)              set -- "$@" "$arg"
  esac
done

while getopts 'oph' opt; do
  case $opt in
    # passthrough the variable to the Makefile for libdash
    o) export optimized_agg_flag=1 ;;
    p) prepare_sudo_install_flag=1 ;;
    *) echo 'Error in command line parsing' >&2
      exit 1
  esac
done
shift "$(( OPTIND - 1 ))"
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
