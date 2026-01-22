[![Build Status](https://travis-ci.com/mgree/libdash.svg?branch=master)](https://travis-ci.com/mgree/libdash)

*libdash* is a fork of the Linux Kernel's `dash` shell that builds a linkable library with extra exposed interfaces. The primary use of libdash is to parse shell scripts, but it could be used for more.

# What are the dependencies?

The C code for dash should build on a wide variety of systems. The library may not build on platforms with esoteric linkers; it's been tested on OS X.

The OCaml code relies on `ctypes-0.11.5` and `ctypes-foreign`; everything else should be in `base`.

# How to build and test it

You should be able to simply run `docker build -t libdash .` to get a runnable environment. Everything will be in `/home/opam/libdash`.

## How to build it locally

Broadly: crib from the `Dockerfile`. More concretely, in the root directory run:

```
./autogen.sh && ./configure && make && sudo make install
```

This should construct an executable `src/dash` and a static library `src/libdash.a`. They will need to be installed globally for things to work well.

Then run:

```
cd ocaml; make && make install
```

This will build the OCaml library and install it in your OPAM repository. There are tests in another directory; they will only build when libdash is actually installed.

```
cd test; make test
```

The tests use `ocaml/round_trip.sh` to ensure that every tester file in `ocaml/tests` round-trips correctly through parsing and pretty printing.

# How to use the parser

The ideal interface to use is `parsecmd_safe` in `parser.c`. Parsing the POSIX shell is a complicated affair: beyond the usual locale issues, aliases affect the lexer, so one must use `setalias` and `unalias` to manage any aliases that ought to exist.

# How work with the parsed nodes

The general AST is described in `nodes.h`. There are some tricky invariants around the precise formatting of control codes; the OCaml code shows some examples of working with the `args` fields in `ocaml/ast.ml`, which converts the C AST to an OCaml AST.
