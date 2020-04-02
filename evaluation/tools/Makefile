.PHONY: all eager-debug

all: eager

eager: eager.c eager_lib.c
	gcc -Wall -O3 eager.c eager_lib.c -o eager

eager-debug: eager.c eager_lib.c
	gcc -Wall -pg -DDEBUG eager.c eager_lib.c -o eager

clean:
	rm eager
