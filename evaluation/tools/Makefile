.PHONY: all eager-debug

all: eager

eager: eager.c
	gcc -Wall -O3 eager.c -o eager

eager-debug: eager.c
	gcc -Wall -pg -DDEBUG eager.c -o eager

clean:
	rm eager
