#!/usr/bin/env bash

# 7zip
sudo apt-get install -y p7zip-full curl wget

# pandoc v.2.2.1
wget https://github.com/jgm/pandoc/releases/download/2.2.1/pandoc-2.2.1-1-$(dpkg --print-architecture).deb
sudo dpkg -i ./pandoc-2.2.1-1-$(dpkg --print-architecture).deb
rm ./pandoc-2.2.1-1-$(dpkg --print-architecture).deb 

# node version 10+ does not need external npm
curl -fsSL https://deb.nodesource.com/setup_10.x | sudo -E bash -
sudo apt-get install -y nodejs
cd  $PASH_TOP/evaluation/scripts/web-index
npm install
cd $PASH_TOP


