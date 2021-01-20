#!/usr/bin/env sh

# TODO: Maybe first check if the repo is accessible via git?

# will install dependencies locally.
PLATFORM=$(uname | tr '[:upper:]' '[:lower:]')
URL='http://get.pash.ndr.md/'
VERSION='latest'
DOWNLOADER='curl'
INSTALL='false'
alias curl='curl -s'

download () {
  command -v curl >/dev/null 2>&1 || { DOWNLOADER='wget'; alias curl='wget -qO- '; }
}

cmd_exists () {
  command -v $1 >/dev/null 2>&1 && echo 'true' || echo 'false';
}

if [ $PLATFORM = 'darwin' ]; then
  echo 'PaSh is not yet well supported on OS X'
fi

while getopts 'p' opt; do
    case $opt in
        p) INSTALL="true" ;;
        *) echo 'Error in command line parsing' >&2
           exit 1
    esac
done
shift "$(( OPTIND - 1 ))"

mkdir pash
cd pash
curl -s $URL | tar xzf - --strip-components=1 -C .
cd scripts
# FIXME: it's unclear we should always run as root
# A first step would be to confirm dependencies are met
# ./install.sh -p
if [[ INSTAll == "true" ]]; then
  echo "ATTEMPT INSTALL"
else
  echo 'Run `./scripts/install.sh`, with `-p` for installing packages if needed'
fi
