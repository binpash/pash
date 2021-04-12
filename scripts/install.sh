#!/usr/bin/env bash

# Ensure that the script fails if something failed
set -e

cd $(dirname $0)
PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}
cd $PASH_TOP

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

## If option -p is set, also run the sudo
if [ "$prepare_sudo_install_flag" -eq 1 ]; then
  ./distro-deps.sh
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