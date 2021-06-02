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
     sudo apt-get update &> $LOG_DIR/apt_update.log
     echo "|-- running apt install..."
     sudo apt-get install -y git libtool m4 automake pkg-config libffi-dev python3 python3-pip curl wamerican-insane bc bsdmainutils \
     samtools bowtie2 vcftools sra-toolkit cutadapt zlib1g-dev default-jre &> $LOG_DIR/apt_install.log
     ;;
   debian*)
     # tested with debian:stable-20210408
     echo "Running preparation sudo apt install:"
     echo "|-- running apt update..."
     apt-get update &> $LOG_DIR/apt_update.log
     echo "|-- running apt install..."
     apt-get install -y git libtool curl sudo procps m4 automake pkg-config curl libffi-dev python3 python3-pip wamerican-insane bc bsdmainutils \
     samtools bowtie2 vcftools sra-toolkit cutadapt zlib1g-dev default-jre &> $LOG_DIR/apt_install.log

     ;;
   fedora*) 
     echo "|-- running dnf install...."
     dnf install -y git gcc python3-pip make automake autoconf libtool hostname \
     bc curl procps samtools bowtie2 vcftools conda zlib-devel java-11-openjdk.x86_64 &> $LOG_DIR/dnf_install.log
     # FIXME sra-toolkit cutadapt for dnf
     ;;
   arch*) 
    echo "Updating mirrors"
    pacman -Sy &> $LOG_DIR/pacman_update.log
     echo "|-- running pacman install...."
    yes | pacman -S git libtool m4 automake pkg-config python-pip libffi make curl autoconf gcc sudo inetutils bc &> $LOG_DIR/pacman_install.log
    # FIXME samtools bowtie2 vcftools conda zlib-devel java-11-openjdk.x86_64 for pacman
    ;;
   *)        echo "unknown distro: '$distro'" ; exit 1 ;;
esac
