#!/bin/bash
IN=../scripts/input/i1M.txt # Change G to M for small input
cat $IN $IN $IN $IN $IN | tr A-Z a-z | grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4'
