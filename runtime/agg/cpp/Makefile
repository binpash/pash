.PHONY: clean
all: wc

CFLAGS = -Wall -std=c++17 -O3 -DNDEBUG -I.
CFLAGS += main.cpp
DEPS = main.cpp main.h

wc: $(DEPS) wc.cpp
	g++ ${CFLAGS} wc.cpp -o wc

clean:
	rm -f wc
