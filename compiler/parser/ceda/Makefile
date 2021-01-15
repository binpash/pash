JSON_C_DIR=/pash/json-c-0.15/install/
SCRIPTS_DIR=/pash/evaluation/

CC=gcc
CFLAGS=-g -O2 -Wall -Werror -I"$(JSON_C_DIR)/include/json-c"
LDFLAGS=-ljson-c -L"$(JSON_C_DIR)/lib" -Wl,-rpath="$(JSON_C_DIR)/lib"


all: json_to_shell2


json_to_shell2: json_to_shell2.o ast2b.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)


tests: json_to_shell2
	@find "$(SCRIPTS_DIR)" | grep '[.]sh$$' | while read f; do sh test_JSON_to_shell2.sh "$$f"; done | tee /tmp/cdash.log
	@echo
	@cat /tmp/cdash.log | grep ':' | cut -d ':' -f 1 | sort | uniq -c


tests-all: json_to_shell2
	@find "/pash/" | grep '[.]sh$$' | while read f; do sh test_JSON_to_shell2.sh "$$f"; done | tee /tmp/cdash.log
	@echo
	@cat /tmp/cdash.log | grep ':' | cut -d ':' -f 1 | sort | uniq -c


clean:
	rm json_to_shell2 *.o
