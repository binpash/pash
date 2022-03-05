#! /usr/bin/env sh
apt-get install software-properties-common -y
add-apt-repository ppa:ubuntu-toolchain-r/test -y
apt-get install g++-10
update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-10 100
update-alternatives --set g++ /usr/bin/g++-10
