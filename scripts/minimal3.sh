#!/bin/bash
IN=./input/i1G.txt # Change G to M for small input
cat $IN | tr A-Z a-z | grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4' | sort
