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
        echo "Running preparation apt install:"
        echo "|-- running apt update..."
        apt-get update &> $LOG_DIR/apt_update.log
        # needed for g++-10
        apt-get install software-properties-common -y &> $LOG_DIR/apt_install.log
        add-apt-repository ppa:ubuntu-toolchain-r/test -y # for g++-10
        echo "|-- running apt install..."
        apt-get install -y git libtool m4 curl automake pkg-config libffi-dev python python3 python3-pip wamerican-insane bc bsdmainutils g++-10 python3-testresources python3-setuptools locales locales-all wget netcat-openbsd &>> $LOG_DIR/apt_install.log
        echo "|-- make g++-10 default..."
        update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-10 100
        update-alternatives --set g++ /usr/bin/g++-10
        ;;
    debian*)
        echo "Running preparation apt install:"
        echo "|-- running apt update..."
        apt-get update &> $LOG_DIR/apt_update.log
        echo "|-- running apt install..."
        apt-get install -y git libtool curl sudo procps m4 automake pkg-config libffi-dev python python3 python3-pip wamerican-insane bc bsdmainutils python3-testresources python3-setuptools netcat-openbsd locales locales-all wget &> $LOG_DIR/apt_install.log
        ;;
    fedora*) 
        echo "|-- running dnf install...."
        dnf install git gcc gcc-c++ python python3-pip make curl automake autoconf libtool hostname bc procps python3-testresources python3-setuptools diffutils python-devel pip python3-setuptools zlib-devel libjpeg-devel nc glibc-langpack-en -y &> $LOG_DIR/dnf_install.log
        ;;
    arch*) 
        echo "Updating mirrors"
        pacman -Sy &> $LOG_DIR/pacman_update.log
        echo "|-- running pacman install...."
        yes | pacman -S git libtool m4 automake curl pkg-config python python-pip libffi make autoconf gcc10 sudo inetutils bc openbsd-netcat
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
