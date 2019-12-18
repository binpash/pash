#!/bin/bash

# Dependencies:
# Python 3.7.3
# The OCaml toplevel, version 4.05.0
# GNU bash, version 5.0.3(1)-release (x86_64-pc-linux-gnu)

python3 --version
ocaml --version
bash --version
# install python 3
pip3 install jsonpickle
sudo apt-get install opam

opam init
eval $(opam env)
opam install dum

git submodule init
git submodule update
git pull

# Setup libdash
cd libdash
./autogen.sh && ./configure && make && sudo make install
cd ocaml
make && sudo make install
dash
cd ../../

# Setup parser
cd parser
make dependencies
make
eval $(opam env)
make clean
cd ../

