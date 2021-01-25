.PHONY: all eager-debug split-debug

all: eager split

CFLAGS=-Wall

ifeq ($(DEBUG),1)
  # We want debug mode.
  CFLAGS += -g -DDEBUG -Og
else
  # We want release mode.
  CFLAGS += -O3 -DNDEBUG
endif

split: split.c
	gcc ${CFLAGS} split.c -o split

split-debug: split.c
	gcc ${CFLAGS} split.c -o split

eager: eager.c eager_lib.c
	gcc ${CFLAGS} eager.c eager_lib.c -o eagers

eager-debug: eager.c eager_lib.c
	gcc ${CFLAGS} -pg eager.c eager_lib.c -o eager

rsplit: r_split.c r_split.h
	gcc ${CFLAGS} r_split.c -o r_split

rmerge: r_merge.c r_split.h
	gcc ${CFLAGS} r_merge.c -o r_merge

rwrap: r_merge.c r_split.h
	gcc ${CFLAGS} r_wrap.c -o r_wrap
clean:
	rm eager split
