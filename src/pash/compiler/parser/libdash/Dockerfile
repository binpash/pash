# start with a reasonable image. Debian 9 stretch is what's on the POSIX testing VM
FROM ocaml/opam2:debian-9

# silence apt
# TODO this still isn't silencing it :(
ENV DEBIAN_FRONTEND=noninteractive

# system support for libdash; libgmp for zarith for lem
RUN sudo apt-get install -y autoconf autotools-dev libtool pkg-config libffi-dev

# because extunix needs camlp4, which isn't ready yet :( 2019-06-24
RUN opam switch 4.07

# make sure we have ocamlfind and ocamlbuild
RUN opam install ocamlfind ocamlbuild

# set up FFI for libdash; num library for lem; extunix for shell syscalls
RUN opam pin add ctypes 0.11.5
RUN opam install ctypes-foreign
RUN opam install extunix

WORKDIR /home/opam

# copy in repo files for libdash to the WORKDIR (should be /home/opam)
# we do this as late as possible so we don't have to redo the slow stuff above
ADD --chown=opam:opam . libdash

# build libdash, expose shared object
RUN cd libdash; ./autogen.sh && ./configure --prefix=/usr --libdir=/usr/lib/x86_64-linux-gnu
RUN cd libdash; make
RUN cd libdash; sudo make install
# build ocaml bindings
RUN cd libdash/ocaml; opam config exec -- make && opam config exec -- make install

# system test
RUN cd libdash/test; opam config exec -- make && opam config exec make test

ENTRYPOINT [ "opam", "config", "exec", "--" ]
CMD [ "bash" ]
