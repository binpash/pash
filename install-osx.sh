
brew install autoconf automake libtool ocaml opam 'python@3.8'

sed 's/libtoolize/glibtoolize/' ./parser/libdash/autogen.sh > autogen-osx.sh
# check if sed command exit correctly before replacing
mv autogen-osx.sh ./parser/libdash/autogen.sh

# Notes:
# if ocamlfind query dash; then ocamlfind remove dash; fi
# ocamlfind: Package `dash' not found

ocamlfind ocamlopt -g -package str,dash,ctypes,ctypes.foreign,atdgen,dum \
        -linkpkg ast_json.ml parse_to_json.ml -o parse_to_json.native
ocamlfind ocamlcp -p a -package str,dash,ctypes,ctypes.foreign,atdgen,dum \
        -linkpkg -i -i ./libdash/ocaml/ast.mli parse_to_json.ml -o parse_to_json.byte
Options -c -o are incompatible with compiling multiple files
make: *** [parse_to_json.byte] Error 2
