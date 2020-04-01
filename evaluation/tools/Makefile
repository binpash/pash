PHONY: all eager-prof

all: eager

eager: eager.c
	gcc -Wall -O3 eager.c -o eager

eager-prof: eager.c
	gcc -Wall -pg eager.c -o eager
