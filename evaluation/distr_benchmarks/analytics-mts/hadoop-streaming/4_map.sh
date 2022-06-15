#!/bin/bash
cat $1 | sed 's/T\(..\):..:../,\1/' | cut -d ',' -f 1,2
