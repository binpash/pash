.PHONY: all eager-debug split-debug

all: eager split

split: split.c
	gcc -Wall -O3 split.c -o split

split-debug: split.c
	gcc -Wall -g -DDEBUG split.c -o split

eager: eager.c eager_lib.c
	gcc -Wall -O3 eager.c eager_lib.c -o eager

eager-debug: eager.c eager_lib.c
	gcc -Wall -pg -DDEBUG eager.c eager_lib.c -o eager

clean:
	rm eager split
