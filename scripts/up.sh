#!/usr/bin/env sh

# clone and setup pash

set -e

# will install dependencies locally.
PLATFORM=$(uname | tr '[:upper:]' '[:lower:]')

if [ "$PLATFORM" = "darwin" ]; then
  echo 'PaSh is not yet well supported on OS X'
  exit 1
fi

set +e
git clone --depth 1 git@github.com:binpash/pash.git
if [ $? -ne 0 ]; then
  echo 'SSH clone failed; attempting HTTPS'
  git clone --depth 1 https://github.com/binpash/pash.git
fi
set -e

if [ $(groups $(whoami) | grep -c "sudo\|root\|admin") -ge 1 ]; then
  # only run this if we are in the sudo group (or it's doomed to fail)
  bash ./pash/scripts/distro-deps.sh
fi
bash ./pash/scripts/setup-pash.sh
