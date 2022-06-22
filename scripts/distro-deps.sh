#!/usr/bin/env bash

# TODO: this should be ran before cloning---otherwise cloning fails.
# It should also set up words etc.

set -e

cd $(dirname $0)
PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}
cd $PASH_TOP

LOG_DIR=$PWD/install_logs
mkdir -p $LOG_DIR

if [[ $(uname) == 'Darwin' ]]; then
  echo 'Currently pash can run only on Linux'
  exit 1
fi


if type lsb_release >/dev/null 2>&1 ; then
   distro=$(lsb_release -i -s)
elif [ -e /etc/os-release ] ; then
   distro=$(awk -F= '$1 == "ID" {print $2}' /etc/os-release)
fi

# convert to lowercase
distro=$(printf '%s\n' "$distro" | LC_ALL=C tr '[:upper:]' '[:lower:]')
# now do different things depending on distro
case "$distro" in
   ubuntu*)  
     echo "Running preparation sudo apt install:"
     echo "|-- running apt update..."
     sudo apt-get update &> $LOG_DIR/apt_update.log
     echo "|-- running apt install..."
     sudo apt-get install -y git libtool m4 curl automake pkg-config libffi-dev python3 python3-pip wamerican-insane bc bsdmainutils &> $LOG_DIR/apt_install.log
     ;;
   debian*)
     # tested with debian:stable-20210408
     echo "Running preparation sudo apt install:"
     echo "|-- running apt update..."
     apt-get update &> $LOG_DIR/apt_update.log
     echo "|-- running apt install..."
     apt-get install -y git libtool curl sudo procps m4 automake pkg-config libffi-dev python3 python3-pip wamerican-insane bc bsdmainutils &> $LOG_DIR/apt_install.log
     ;;
   fedora*) 
     echo "|-- running dnf install...."
     dnf install git gcc python3-pip make curl automake autoconf libtool hostname bc procps -y  &> $LOG_DIR/dnf_install.log
     ;;
   arch*) 
    echo "Updating mirrors"
    pacman -Sy &> $LOG_DIR/pacman_update.log
     echo "|-- running pacman install...."
    yes | pacman -S git libtool m4 automake curl pkg-config python-pip libffi make autoconf gcc sudo inetutils bc
    ;;
   *)        echo "unknown distro: '$distro'" ; exit 1 ;;
esac
