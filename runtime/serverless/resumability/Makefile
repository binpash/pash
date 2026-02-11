all: out-gate-stateful out-gate-stateless in-gate-stateful in-gate-stateless in-gate-stateful-readline
.PHONY: all clean

CFLAGS=-Wall

out-gate-stateful: out_gate_stateful.c r_split.h
	gcc ${CFLAGS} out_gate_stateful.c -o out_gate_stateful

out-gate-stateless: out_gate_stateless.c r_split.h
	gcc ${CFLAGS} out_gate_stateless.c -o out_gate_stateless

in-gate-stateless: in_gate_stateless.c r_split.h
	gcc ${CFLAGS} in_gate_stateless.c -o in_gate_stateless

in-gate-stateful: in_gate_stateful.c r_split.h
	gcc ${CFLAGS} in_gate_stateful.c -o in_gate_stateful

in-gate-stateful-readline: in_gate_stateful_readline.c r_split.h
	gcc ${CFLAGS} in_gate_stateful_readline.c -o in_gate_stateful_readline

clean:
	rm -f r_split r_merge r_wrap r_unwrap out_gate in_gate out_gate_stateful in_gate_stateful out-gate-manual
