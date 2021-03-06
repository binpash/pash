## Instructions

The parser in this directory uses the mgree/libdash posix compliant
parser and outputs the AST in JSON format using atdgen.

In order to install, one has to execute `make opam-dependencies && make libdash && make`. The first command makes libdash and sets up all the ocaml dependencies.

Alternatively, one can install libdash (as explained in its README)
and then run `make` here.

To run the parser, one can run:

```sh
./parse_to_json.native <sh file>
```

The following process was followed to make the parser output the ast to json.

* specify the AST definition in `ast_atd.atd`.
* ran `atdgen -j -j-std ast_atd.atd` to produce `ast_atd_j.ml`.
* copy it to `ast_json.ml` (removing the `char` definition).
* also make a small adjustment in the `ast.ml`, `ast.mli` files in `libdash/ocaml`. (?nv)

This procedure is not automated, and in case the AST definition in
libdash changes, this process has to be done again.
