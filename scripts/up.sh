#!/usr/bin/env sh

# clone and setup pash

set -e

# will install dependencies locally.
PLATFORM=$(uname | tr '[:upper:]' '[:lower:]')

if [ "$PLATFORM" = "darwin" ]; then
  echo 'PaSh is not yet well supported on OS X'
  exit 1
fi

## Download the latest PaSh release tarball to avoid using git
curl -s https://api.github.com/repos/binpash/pash/releases/latest | 
  grep "tarball_url" | 
  cut -d : -f 2,3 | 
  tr -d \" | 
  tr -d , | 
  wget -i - -O pash.tar.gz
## Find the name of the top directory
pash_dir_name=`tar -tzf pash.tar.gz | head -1 | cut -f1 -d"/"`
tar -xzf pash.tar.gz
mv "$pash_dir_name" pash
# Alternative:
# git clone https://github.com/binpash/pash.git

if [ $(groups $(whoami) | grep -c "sudo\|root\|admin") -ge 1 ]; then
  # only run this if we are in the sudo group (or it's doomed to fail)
  bash ./pash/scripts/distro-deps.sh
fi
export PASH_TOP="$PWD/pash"
bash ./pash/scripts/setup-pash.sh
