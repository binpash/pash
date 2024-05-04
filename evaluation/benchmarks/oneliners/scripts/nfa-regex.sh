#!/bin/bash
# Match complex regular-expression over input

cd "$(dirname "$0")" || exit 1

SIZE=200M

IN=input/$SIZE.txt

cat $IN | tr A-Z a-z | grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4'
