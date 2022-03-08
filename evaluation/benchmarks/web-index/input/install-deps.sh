#!/usr/bin/env bash

# 7zip
pkgs='p7zip-full curl wget nodejs' 
if ! dpkg -s $pkgs >/dev/null 2>&1 ; then
  sudo apt-get install $pkgs -y
  echo 'Packages Installed'
fi

if ! dpkg -s pandoc > /dev/null 2>&1 ; then
  # pandoc v.2.2.1
  wget https://github.com/jgm/pandoc/releases/download/2.2.1/pandoc-2.2.1-1-$(dpkg --print-architecture).deb
  sudo dpkg -i ./pandoc-2.2.1-1-$(dpkg --print-architecture).deb
  rm ./pandoc-2.2.1-1-$(dpkg --print-architecture).deb 
fi

if [ ! -d node_modules ]; then
  # node version 10+ does not need external npm
  curl -fsSL https://deb.nodesource.com/setup_10.x | sudo -E bash -
  npm install
fi
