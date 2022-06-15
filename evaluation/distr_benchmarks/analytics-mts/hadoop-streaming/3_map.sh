#!/bin/bash
cat $IN | sed 's/T\(..\):..:../,\1/' |  cut -d ',' -f 1,2,4
