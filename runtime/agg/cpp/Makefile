.PHONY: clean
all: wc uniq-c

CFLAGS = -Wall -std=c++17 -O3 -DNDEBUG -I.
CFLAGS += main.cpp
DEPS = main.cpp main.h

wc: $(DEPS) wc.cpp
	g++ ${CFLAGS} wc.cpp -o wc
uniq-c: $(DEPS) uniq-c.cpp
	g++ ${CFLAGS} uniq-c.cpp -o uniq-c

clean:
	rm -f wc uniq-c
