PARSER=parse_to_json
JSON_TO_SHELL=json_to_shell
AST_I=./libdash/ocaml/ast.mli

.PHONY : all test dependencies libdash libdash-ocaml clean

all : $(PARSER).native $(PARSER).byte $(JSON_TO_SHELL).native $(JSON_TO_SHELL).byte

test : $(PARSER).native $(PARSER).byte $(wildcard tests/*)
	@echo "TESTING test.native"
	@for f in tests/*; do \
		./round_trip.sh ./$(PARSER).native $$f 2>test.err; \
	done
	@echo "TESTING test.byte"
	@for f in tests/*; do \
		./round_trip.sh ./$(PARSER).byte $$f 2>test.err; \
	done

## TODO: Unify the two rules below
$(PARSER).native : ast_json.ml $(PARSER).ml 
	ocamlfind ocamlopt -g -package str,dash,ctypes,ctypes.foreign,atdgen,dum \
	-linkpkg $^ -o $(PARSER).native

$(PARSER).byte : $(PARSER).ml 
	ocamlfind ocamlcp -p a -package str,dash,ctypes,ctypes.foreign,atdgen,dum \
	-linkpkg -i -i $(AST_I) $^ -o $(PARSER).byte

$(JSON_TO_SHELL).native : ast_json.ml $(JSON_TO_SHELL).ml 
	ocamlfind ocamlopt -g -package str,dash,ctypes,ctypes.foreign,atdgen,dum \
	-linkpkg $^ -o $(JSON_TO_SHELL).native

$(JSON_TO_SHELL).byte : $(JSON_TO_SHELL).ml 
	ocamlfind ocamlcp -p a -package str,dash,ctypes,ctypes.foreign,atdgen,dum \
	-linkpkg -i -i $(AST_I) $^ -o $(JSON_TO_SHELL).byte

clean :
	rm -f *.o *.cmo *.cmi *.cmx $(PARSER).native $(PARSER).byte \
	$(JSON_TO_SHELL).native $(JSON_TO_SHELL).byte test.err

libdash:
	cd libdash && ./autogen.sh || ./autogen.sh && ./configure && make && sudo make install

libdash-ocaml:
	cd libdash/ocaml && make && make install

opam-dependencies: 
	opam install -y ctypes.0.11.5
	opam install -y ctypes-foreign
	opam install -y atdgen
	opam install -y dum
