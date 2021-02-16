.PHONY: all eager-debug split-debug

all: eager split r-merge r-wrap r-split r-unwrap

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
  CXXFLAGS += -g -ggdb
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
clean:
	rm eager split
