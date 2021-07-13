.PHONY: clean
all: grep grep-c

CFLAGS = -Wall -std=c++17 -O3 -DNDEBUG -I.
CFLAGS += main.cpp
DEPS = main.cpp main.h

grep: $(DEPS) grep.cpp
	g++ ${CFLAGS} grep.cpp -o grep
grep-c: $(DEPS) grep-c.cpp
	g++ ${CFLAGS} grep-c.cpp -o grep-c

clean:
	rm -f grep grep-c
