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

echo "Running preparation sudo apt install:"
echo "|-- running apt update..."
sudo apt-get update &> $LOG_DIR/apt_update.log
echo "|-- running apt install..."
sudo apt-get install -y git libtool m4 automake pkg-config libffi-dev python3 python3-pip wamerican-insane bc bsdmainutils &> $LOG_DIR/apt_install.log
