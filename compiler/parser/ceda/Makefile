JSON_C_DIR=/pash/json-c-0.15/install/
SCRIPTS_DIR=/pash/

CC=gcc
CFLAGS=-g -O2 -Wall -Werror -I"$(JSON_C_DIR)/include/json-c" -I"../libdash/src"
LDFLAGS=-ljson-c -L"$(JSON_C_DIR)/lib" -Wl,-rpath="$(JSON_C_DIR)/lib"


all: parse_to_json2 json_to_shell2 prettyprint_json


parse_to_json2: parse_to_json2.o ast2a.o ast2json.o dash2.o CharList.o Stack.o arg_char.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) -ldash

json_to_shell2: json_to_shell2.o ast2b.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

prettyprint_json: prettyprint_json.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)


tests: parse_to_json2 json_to_shell2
	@find "$(SCRIPTS_DIR)" | grep '[.]sh$$' | while read f; do sh test_rt.sh "$$f"; done | tee /tmp/cdash.log
	@echo
	@cat /tmp/cdash.log | egrep '^[A-Z]+:' | cut -d ':' -f 1 | sort | uniq -c

testsA: parse_to_json2
	@find "$(SCRIPTS_DIR)" | grep '[.]sh$$' | while read f; do sh test_parse_to_JSON2.sh "$$f"; done | tee /tmp/cdashA.log
	@echo
	@cat /tmp/cdashA.log | egrep '^[A-Z]+:' | cut -d ':' -f 1 | sort | uniq -c

tests: test_rt.sh parse_to_json2 json_to_shell2
	@find "$(SCRIPTS_DIR)" | grep '[.]sh$$' | while read f; do sh test_rt.sh "$$f"; done | tee /tmp/cdash_rt.log
	@echo
	@cat /tmp/cdash_rt.log | egrep '^[A-Z]+:' | cut -d ':' -f 1 | sort | uniq -c


testsA_py: test_rt_py.sh
	@find "$(SCRIPTS_DIR)" | grep '[.]sh$$' | while read f; do sh test_rt_py.sh "$$f"; done | tee /tmp/cdashA_py.log
	@echo
	@cat /tmp/cdashA_py.log | egrep '^[A-Z]+:' | cut -d ':' -f 1 | sort | uniq -c


clean:
	rm parse_to_json2 json_to_shell2 *.o
