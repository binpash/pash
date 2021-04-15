#!/bin/bash

custom_sort() {
    sort $@
}

custom_tr() {
    tr A-Z a-z
}

export -f custom_tr

FILES="../evaluation/tests/input/1M.txt ../evaluation/tests/input/1M.txt"

cat $FILES | custom_tr | custom_sort
