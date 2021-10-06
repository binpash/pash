#!/usr/bin/env bash

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
     sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y # for g++-10
     sudo apt-get update &> $LOG_DIR/apt_update.log
     echo "|-- running apt install..."
     sudo apt-get install -y git libtool m4 curl automake pkg-config libffi-dev python3 python3-pip wamerican-insane bc bsdmainutils g++-10 python3-testresources &> $LOG_DIR/apt_install.log
     echo "|-- make g++-10 default..."
     sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-10 100
     sudo update-alternatives --set g++ /usr/bin/g++-10
     ;;
   debian*)
     # tested with debian:stable-20210408
     echo "Running preparation sudo apt install:"
     echo "|-- running apt update..."
     add-apt-repository ppa:ubuntu-toolchain-r/test -y # for g++-10
     apt-get update &> $LOG_DIR/apt_update.log
     echo "|-- running apt install..."
     apt-get install -y git libtool curl sudo procps m4 automake pkg-config libffi-dev python3 python3-pip wamerican-insane bc bsdmainutils g++-10 python3-testresources &> $LOG_DIR/apt_install.log
     echo "|-- make g++-10 default..."
     sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-10 100
     sudo update-alternatives --set g++ /usr/bin/g++-10
     ;;
   fedora*) 
     echo "|-- running dnf install...."
     dnf install git gcc gcc-c++ python3-pip make curl automake autoconf libtool hostname bc procps python3-testresources -y &> $LOG_DIR/dnf_install.log
     ;;
   arch*) 
     echo "Updating mirrors"
     pacman -Sy &> $LOG_DIR/pacman_update.log
     echo "|-- running pacman install...."
     yes | pacman -S git libtool m4 automake curl pkg-config python-pip libffi make autoconf gcc10 sudo inetutils bc
     ;;
   freebsd*)
     echo "Updating mirros"
     pkg update &> $LOG_DIR/pkg_update.log
     echo "|-- running pkg install...."
     # TODO add python3-testresources dep
     yes | pkg install libtool m4 automake curl libffi py38-pip autoconf gcc gsed gmake
     ;;
   *)        echo "unknown distro: '$distro'" ; exit 1 ;;
esac
