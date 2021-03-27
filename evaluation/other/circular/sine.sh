#!/bin/bash
F="temp.txt"
[ -f $F ] && (rm $F && echo 1 >$F )
tail -f temp.txt | while read n; do echo "1+s(3*$n)" | bc -l; sleep 1; done | tee -a temp.txt
