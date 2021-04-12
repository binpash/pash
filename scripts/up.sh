#!/usr/bin/env sh

set -e

# will install dependencies locally.
PLATFORM=$(uname | tr '[:upper:]' '[:lower:]')
URL='https://github.com/andromeda/pash/archive/refs/heads/main.zip'
VERSION='latest'
DL=$(command -v curl >/dev/null 2>&1 && echo curl || echo 'wget -qO-')

cmd_exists () {
  command -v $1 >/dev/null 2>&1 && echo 'true' || echo 'false';
}

if [[ "$PLATFORM" == "darwin" ]]; then
  echo 'PaSh is not yet well supported on OS X'
  exit 1
fi

git clone git@github.com:andromeda/pash.git
cd pash/scripts
./distro-deps.sh
./setup-pash.sh

