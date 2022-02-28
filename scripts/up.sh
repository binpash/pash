#!/usr/bin/env sh

# clone and setup pash
# N.b. This is a .sh script

set -e

# will install dependencies locally.
PLATFORM=$(uname | tr '[:upper:]' '[:lower:]')
URL='https://github.com/binpash/pash/archive/refs/heads/main.zip'
VERSION='latest'
DL=$(command -v curl >/dev/null 2>&1 && echo curl || echo 'wget -qO-')

cmd_exists () {
  command -v $1 >/dev/null 2>&1 && echo 'true' || echo 'false';
}

if [ "$PLATFORM" = "darwin" ]; then
  echo 'PaSh is not yet well supported on OS X'
  exit 1
fi

set +e
git clone git@github.com:binpash/pash.git
if [ $? -ne 0 ]; then
  echo 'SSH clone failed; attempting HTTPS'
  git clone https://github.com/andromeda/pash.git
fi
set -e

cd pash/scripts
# git checkout s3 # FIXME only for testing while PR is up

if [ $(groups $(whoami) | grep -c "sudo\|root\|admin") -ge 1 ]; then
  # only run this if we are in the sudo group (or it's doomed to fail)
  bash distro-deps.sh
fi
bash setup-pash.sh
