#!/usr/bin/env bash

set -e

cd $(dirname $0)
PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}
optimized_agg_flag=${optimized_agg_flag:-0}
show_deps=${show_deps:-0}
cd $PASH_TOP

LOG_DIR=install_logs
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
        pkgs="git libtool m4 curl automake pkg-config libffi-dev python python3 python3-pip wamerican-insane bc bsdmainutils python3-testresources python3-setuptools locales locales-all wget netcat-openbsd" 
        if [[ "$show_deps" == 1 ]]; then
            echo "$pkgs"
            exit 0
        fi      
        echo "Running preparation apt install:"
        echo "|-- running apt update..."
        apt-get update &> $LOG_DIR/apt_update.log
        echo "|-- running apt install..."
        apt-get install -y $pkgs &>> $LOG_DIR/apt_install.log
        if [[ "$optimized_agg_flag" == 1 ]];  then
            echo "|-- installing g++-10..."
            apt-get install software-properties-common -y &> $LOG_DIR/apt_install.log
            add-apt-repository ppa:ubuntu-toolchain-r/test -y  &> $LOG_DIR/apt_install.log
            apt-get install g++-10  &> $LOG_DIR/apt_install.log
            update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-10 100 &> $LOG_DIR/apt_install.log
            update-alternatives --set g++ /usr/bin/g++-10 &> $LOG_DIR/apt_install.log
        fi
        ;;
    debian*)
        pkgs="git libtool curl sudo procps m4 automake pkg-config libffi-dev python python3 python3-pip wamerican-insane bc bsdmainutils python3-testresources python3-setuptools netcat-openbsd locales locales-all wget"
        if [[ "$show_deps" == 1 ]]; then
            echo "$pkgs"
            exit 0
        fi
        echo "Running preparation apt install:"
        echo "|-- running apt update..."
        apt-get update &> $LOG_DIR/apt_update.log
        echo "|-- running apt install..."
        apt-get install -y $pkgs &> $LOG_DIR/apt_install.log
        ;;
    fedora*) 
        pkgs="git gcc gcc-c++ python python3-pip make curl automake autoconf libtool hostname bc procps python3-testresources python3-setuptools diffutils python-devel pip python3-setuptools zlib-devel libjpeg-devel nc glibc-langpack-en"
        if [[ "$show_deps" == 1 ]]; then
            echo "$pkgs"
            exit 0
        fi
        echo "|-- running dnf install...."
        dnf install -y $pkgs &> $LOG_DIR/dnf_install.log
        ;;
    arch*) 
        pkgs="git libtool m4 automake curl pkg-config python python-pip libffi make autoconf gcc10 sudo inetutils bc openbsd-netcat"
        if [[ "$show_deps" == 1 ]]; then
            echo "$pkgs"
            exit 0
        fi
        echo "Updating mirrors"
        pacman -Sy &> $LOG_DIR/pacman_update.log
        echo "|-- running pacman install...."
        yes | pacman -S $pkgs &> $LOG_DIR/pacman_install.log
        ;;
    freebsd*)
        pkgs="libtool m4 automake curl libffi py38-pip autoconf gcc gsed gmake"
        if [[ "$show_deps" == 1 ]]; then
            echo "$pkgs"
            exit 0
        fi
        echo "Updating mirrors"
        pkg update &> $LOG_DIR/pkg_update.log
        echo "|-- running pkg install...."
        # TODO add python3-testresources dep
        yes | pkg install $pkgs
        ;;
    *)        echo "unknown distro: '$distro'" ; exit 1 ;;
esac
