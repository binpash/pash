#!/bin/sh

# MIT License

# Copyright (c) 2017 Stockholm University

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

install_dir=/opt/nmrpipe

set -e # fail as soon there is an error

# check that we are running a 64 bit system

num_bits=`getconf LONG_BIT`
if $num_bits != "64"; then
  echo This is not a 64 bit OS
  exit 1
fi

download_dir=`mktemp -d /tmp/tmp_download.XXX`


URL_list="https://www.ibbr.umd.edu/nmrpipe/install.com \
          https://www.ibbr.umd.edu/nmrpipe/binval.com \
          https://www.ibbr.umd.edu/nmrpipe/NMRPipeX.tZ \
          https://www.ibbr.umd.edu/nmrpipe/s.tZ \
          https://www.ibbr.umd.edu/nmrpipe/dyn.tZ \
          https://www.ibbr.umd.edu/nmrpipe/talos.tZ \
	  http://spin.niddk.nih.gov/bax/software/smile/plugin.smile.tZ"

for i in $URL_list; do
    echo downloading $i
    wget "--directory-prefix=$download_dir" --no-directories $i
done

sudo mkdir "$install_dir"
sudo find $download_dir -maxdepth 1 -mindepth 1 -type f -name '*.com' -exec cp {} "$install_dir" \;

sudo find "$install_dir" -maxdepth 1 -mindepth 1 -type f -name '*.com' -exec chmod 755 {} \;

sudo dpkg --add-architecture i386

sudo apt-get update

echo ttf-mscorefonts-installer msttcorefonts/accepted-mscorefonts-eula select true | sudo debconf-set-selections

sudo apt-get install -y csh
sudo apt-get install -y default-jdk
sudo apt-get install -y default-jre
sudo apt-get install -y libc6:i386
sudo apt-get install -y libstdc++6:i386
sudo apt-get install -y libx11-6:i386
sudo apt-get install -y libxext6:i386
sudo apt-get install -y msttcorefonts

sudo sh -c "cd $install_dir ; ./install.com +src $download_dir option +type linux212_64"
