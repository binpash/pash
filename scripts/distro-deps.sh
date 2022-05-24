#!/usr/bin/env bash

set -e
cd $(dirname $0)
PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}
. "$PASH_TOP/scripts/utils.sh"

if [[ $(uname) == 'Darwin' ]]; then
    echo 'Currently pash can run only on Linux'
    exit 1
fi

read_cmd_args $@
cd $PASH_TOP

LOG_DIR=install_logs
mkdir -p $LOG_DIR

# if we aren't running in docker, use sudo to install packages
if [ ! -f /.dockerenv ]; then
    export SUDO="sudo"
fi  

if type lsb_release >/dev/null 2>&1 ; then
    distro=$(lsb_release -i -s)
elif [ -e /etc/os-release ] ; then
    distro=$(awk -F= '$1 == "ID" {print $2}' /etc/os-release)
fi

# convert to lowercase
distro=$(printf '%s\n' "$distro" | LC_ALL=C tr '[:upper:]' '[:lower:]')
# compile the list of the shared required packages
pkgs="automake bc curl gcc git graphviz libtool m4 python sudo wget"
# now do different things depending on distro
case "$distro" in
    ubuntu*)  
        pkgs="$pkgs bsdmainutils libffi-dev locales locales-all netcat-openbsd pkg-config python3 python3-pip python3-setuptools python3-testresources wamerican-insane"
        if [[ "$show_deps" == 1 ]]; then
            echo "$pkgs" | sort
            exit 0
        fi      
        echo "Running preparation apt install:"
        echo "|-- running apt update..."
        $SUDO apt-get update &> $LOG_DIR/apt_update.log
        echo "|-- running apt install..."
        $SUDO apt-get install -y $pkgs &>> $LOG_DIR/apt_install.log
        if [[ "$optimized_agg_flag" == 1 ]];  then
            echo "|-- installing g++-10..."
            $SUDO apt-get install software-properties-common -y &> $LOG_DIR/apt_install.log
            $SUDO add-apt-repository ppa:ubuntu-toolchain-r/test -y  &> $LOG_DIR/apt_install.log
            $SUDO apt-get install g++-10 -y &> $LOG_DIR/apt_install.log
            $SUDO update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-10 100 &> $LOG_DIR/apt_install.log
            $SUDO update-alternatives --set g++ /usr/bin/g++-10 &> $LOG_DIR/apt_install.log
        fi
        ;;
    debian*)
        pkgs="$pkgs bsdmainutils libffi-dev locales locales-all netcat-openbsd pkg-config procps python3 python3-pip python3-setuptools python3-testresources wamerican-insane"
        if [[ "$show_deps" == 1 ]]; then
            echo "$pkgs" | sort
            exit 0
        fi
        echo "Running preparation apt install:"
        echo "|-- running apt update..."
        $SUDO apt-get update &> $LOG_DIR/apt_update.log
        echo "|-- running apt install..."
        $SUDO apt-get install -y $pkgs &> $LOG_DIR/apt_install.log
        ;;
    fedora*) 
        pkgs="$pkgs autoconf diffutils gcc-c++ glibc-langpack-en hostname libjpeg-devel make nc pip procps python-devel python3-pip python3-setuptools python3-setuptools python3-testresources zlib-devel"
        if [[ "$show_deps" == 1 ]]; then
            echo "$pkgs" | sort
            exit 0
        fi
        echo "|-- running dnf install...."
        $SUDO dnf install -y $pkgs &> $LOG_DIR/dnf_install.log
        ;;
    arch*) 
        pkgs="$pkgs autoconf inetutils libffi make openbsd-netcat pkg-config python-pip"
        if [[ "$show_deps" == 1 ]]; then
            echo "$pkgs" | sort
            exit 0
        fi
        echo "Updating mirrors"
        $SUDO pacman -Sy &> $LOG_DIR/pacman_update.log
        echo "|-- running pacman install...."
        yes | $SUDO pacman -S $pkgs &> $LOG_DIR/pacman_install.log
        ;;
    freebsd*)
        pkgs="$pkgs autoconf gmake gsed libffi py38-pip"
        if [[ "$show_deps" == 1 ]]; then
            echo "$pkgs" | sort
            exit 0
        fi
        echo "Updating mirrors"
        $SUDO pkg update &> $LOG_DIR/pkg_update.log
        echo "|-- running pkg install...."
        # TODO add python3-testresources dep
        yes | $SUDO pkg install $pkgs
        ;;
    *)        echo "unknown distro: '$distro'" ; exit 1 ;;
esac
