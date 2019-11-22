#!/bin/bash

# A script that showcases truly async pipes (via fs)
# Note to self: remember | { lambda }

fz () { sleep $0; echo "1-"$0; }

export -f fz

: > f1

tail -f ./f1 | xargs -n 1  bash -c 'fz "$@"' &

# {seq 5; echo 'yay!' >&2 ; } > ./f1
seq 5 > ./f1
echo 'yay!'
