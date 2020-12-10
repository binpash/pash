#!/bin/bash

# TODO: Maybe first check if the repo is accessible via git?

# will install dependencies locally.
PLATFORM=$(uname | tr '[:upper:]' '[:lower:]')
URL='http://get.pash.ndr.md/'
VERSION='latest'
DOWNLOADER='curl'
alias curl='curl -s'

download () {
    command -v curl >/dev/null 2>&1 || 
        { DOWNLOADER='wget'; alias curl='wget -qO- '; }
}

cmd_exists () {
    command -v $1 >/dev/null 2>&1 && 
        echo 'true' ||
        echo 'false';
}

if [ $PLATFORM = 'darwin' ]; then
  echo 'PaSh is not yet well supported on OS X'
fi

curl -s $URL | tar xzf - --strip-components=1 -C .
cd pash/scripts
./install.sh
