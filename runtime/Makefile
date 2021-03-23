.PHONY: all eager-debug split-debug

all: eager split r-merge r-wrap r-split r-unwrap dgsh-tee

CFLAGS=-Wall

ifeq ($(DEBUG),1)
  # We want debug mode.
  CFLAGS += -g -DDEBUG -Og -ggdb -Wno-attributes
else
  # We want release mode.
  CFLAGS += -O3 -DNDEBUG
endif

ifeq ($(PERF),1)
  # We want debug mode.
  CFLAGS += -g -ggdb
endif

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S), Darwin)
DGSH_ASSEMBLY_FILE=dgsh-macho.s
else
DGSH_ASSEMBLY_FILE=dgsh-elf.s
endif


split: split.c
	gcc ${CFLAGS} split.c -o split

split-debug: split.c
	gcc ${CFLAGS} split.c -o split

eager: eager.c eager_lib.c
	gcc ${CFLAGS} eager.c eager_lib.c -o eager

eager-debug: eager.c eager_lib.c
	gcc ${CFLAGS} -pg eager.c eager_lib.c -o eager

r-split: r_split.c r_split.h
	gcc ${CFLAGS} r_split.c -o r_split

r-merge: r_merge.c r_split.h
	gcc ${CFLAGS} r_merge.c -o r_merge

r-wrap: r_wrap.c r_split.h
	gcc ${CFLAGS} r_wrap.c -o r_wrap

r-unwrap: r_unwrap.c r_split.h
	gcc ${CFLAGS} r_unwrap.c -o r_unwrap


libdgsh_a_SOURCES = negotiate.c $(DGSH_ASSEMBLY_FILE)
dgsh_tee_SOURCES = dgsh-tee.c

## For some reason needs DDEBUG flag to work properly with testing script
## TODO: figure out how to make it work without DDEBUG
dgsh-tee:
	git clone https://github.com/dspinellis/dgsh
	cd dgsh/core-tools/src; \
	git checkout 9e695a59c3d6ce7de14fdb6d9f9613454cf7281b;\
	gcc ${CFLAGS} $(dgsh_tee_SOURCES) $(libdgsh_a_SOURCES) -o $@; \
	cp dgsh-tee ../../../dgsh-tee;
	rm -rf dgsh/

clean:
	rm -f eager split r_split r_wrap r_unwrap dgsh-tee
	rm -rf dgsh