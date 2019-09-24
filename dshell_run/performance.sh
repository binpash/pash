#!/bin/bash

COUNTER=0
while [ $COUNTER -lt 10000 ]; do
    cat /home/cvetkovic/sdsh/scripts/input.txt | wc -w
    let COUNTER=COUNTER+1
done
