.PHONY: clean
all: grep grep-c

CFLAGS = -Wall -O3 -DNDEBUG -I.
CFLAGS += main.cpp
DEPS = main.cpp main.h

grep: $(DEPS) grep.cpp
	gcc ${CFLAGS} grep.cpp -o grep
grep-c: $(DEPS) grep-c.cpp
	gcc ${CFLAGS} grep-c.cpp -o grep-c

clean:
	rm -f grep grep-c
