#!/bin/bash

custom_sort() {
    sort $@
}

FILES="../evaluation/tests/input/1M.txt ../evaluation/tests/input/1M.txt"

cat $FILES | tr A-Z a-z | custom_sort
