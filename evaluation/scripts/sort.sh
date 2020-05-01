#!/bin/bash

IN=../evaluation/scripts/input/10G.txt

cat $IN | tr A-Z a-z | sort --parallel=$1
