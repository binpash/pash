#!/bin/bash

custom_sort() {
    sort $@
}

FILES="../evaluation/scripts/input/1M.txt ../evaluation/scripts/input/1M.txt"

cat $FILES | tr A-Z a-z | custom_sort
