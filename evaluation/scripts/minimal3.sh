#!/bin/bash
IN=../scripts/input/i10M.txt # Change G to M for small input
cat $IN $IN | tr A-Z a-z | grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4'
