#!/bin/bash

custom_sort() {
    sort $@
}

FILE=../evaluation/scripts/input/1M.txt

cat $FILE | tr A-Z a-z | custom_sort
