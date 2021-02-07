JSON_C_DIR=/pash/json-c-0.15/install/
SCRIPTS_DIR=/pash/

CC=gcc
CFLAGS=-fPIC -g -O2 -Wall -Werror -I"$(JSON_C_DIR)/include/json-c" -I"../libdash/src"
LDFLAGS=-ljson-c -L"$(JSON_C_DIR)/lib" -Wl,-rpath="$(JSON_C_DIR)/lib"


all: parse_to_json2 json_to_shell2 prettyprint_json libdash2.so


# C binaries
parse_to_json2: parse_to_json2.o ast2a.o ast2json.o dash2.o CharList.o Stack.o arg_char.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) -ldash

json_to_shell2: json_to_shell2.o ast2b.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

prettyprint_json: prettyprint_json.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

# Shared objects
libdash2.so: dash2.o CharList.o
	$(CC) -shared $^ -o $@ $(LDFLAGS) -ldash

# Tests for C binaries
testsA: parse_to_json2
	@find "$(SCRIPTS_DIR)" | grep '[.]sh$$' | while read f; do sh test_parse_to_JSON2.sh "$$f"; done | tee /tmp/cdashA.log
	@echo
	@cat /tmp/cdashA.log | egrep '^[A-Z0-9_]+:' | cut -d ':' -f 1 | sort | uniq -c

testsB: json_to_shell2
	@find "$(SCRIPTS_DIR)" | grep '[.]sh$$' | while read f; do sh test_JSON_to_shell2.sh "$$f"; done | tee /tmp/cdashB.log
	@echo
	@cat /tmp/cdashB.log | egrep '^[A-Z0-9_]+:' | cut -d ':' -f 1 | sort | uniq -c

testsRT: test_rt.sh parse_to_json2 json_to_shell2
	@find "$(SCRIPTS_DIR)" | grep '[.]sh$$' | while read f; do sh test_rt.sh "$$f"; done | tee /tmp/cdash_rt.log
	@echo
	@cat /tmp/cdash_rt.log | egrep '^[A-Z0-9_]+:' | cut -d ':' -f 1 | sort | uniq -c

# Tests for Python scripts
testsRT_py: test_rt_py.sh
	@find "$(SCRIPTS_DIR)" | grep '[.]sh$$' | while read f; do sh test_rt_py.sh "$$f"; done | tee /tmp/cdash_rt_py.log
	@echo
	@cat /tmp/cdash_rt_py.log | egrep '^[A-Z0-9_]+:' | cut -d ':' -f 1 | sort | uniq -c

#testsB_py: test_ast2shell_py.sh
#	@find "$(SCRIPTS_DIR)" | grep '[.]sh$$' | while read f; do sh test_ast2shell_py.sh "$$f"; done | tee /tmp/cdashB_py.log
#	@echo
#	@cat /tmp/cdashB_py.log | egrep '^[A-Z]+:' | cut -d ':' -f 1 | sort | uniq -c

# Clean
clean:
	rm parse_to_json2 json_to_shell2 *.o *.so
