#!/usr/bin/env sh

# clone and setup pash
# N.b. This is a .sh script

set -e

## will install dependencies locally.
#PLATFORM=$(uname | tr '[:upper:]' '[:lower:]')
#URL='https://github.com/andromeda/pash/archive/refs/heads/main.zip'
#VERSION='latest'
#DL=$(command -v curl >/dev/null 2>&1 && echo curl || echo 'wget -qO-')
#
#cmd_exists () {
#  command -v $1 >/dev/null 2>&1 && echo 'true' || echo 'false';
#}
#
#if [ "$PLATFORM" = "darwin" ]; then
#  echo 'PaSh is not yet well supported on OS X'
#  exit 1
#fi
#
#git clone git@github.com:andromeda/pash.git
cd pash/scripts
# git checkout s3 # FIXME only for testing while PR is up
bash distro-deps.sh
bash setup-pash.sh
