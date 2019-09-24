#!/bin/bash

cat ./input/i1G.txt | tr A-Z a-z | grep '\(.\).*\1\(.\).*\2\(.\).*\3\(.\).*\4' > /dev/null
