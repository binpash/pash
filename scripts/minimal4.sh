#!/bin/bash
IN=../scripts/input/i1G.txt # Change G to M for small input
cat $IN $IN | tr A-Z a-z | tr a-z A-Z
