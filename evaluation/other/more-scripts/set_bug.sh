#!/bin/bash

f() {
    echo "f: $@"
}

set -- a b c
echo "$@"
f
echo "$@"

