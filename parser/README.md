## Instructions

The parser in this directory uses the mgree/libdash posix compliant
parser and outputs the AST in JSON format using atdgen.

In order for this to work, one has to install the libdash (as
explained in its README) and then run `make` here.

To run the parser, one can run:

```sh
./parse_to_json.native <sh file>
```

The following process was followed to make the parser output the ast to json.

1. I specified the AST definition in `ast_atd.atd`.

2. I ran `atdgen -j ast_atd.atd` which produced the `ast_atd_j.ml` file.

3. I then copied the part of the file without the type definitions to `ast_json.ml`.

4. I also had to make a small adjustment in the `ast.ml`, `ast.mli` files in `libdash/ocaml`.

This procedure is not automated, and in case the AST definition in
libdash changes, this process has to be done again. However, this
should be enough for our purpose now.
