#!/bin/bash
# https://unix.stackexchange.com/questions/193441/how-can-i-implement-a-circular-flow-of-data-among-interconnected-commands
F="temp.txt"
[ -f $F ] && (rm $F && echo 1 >$F )
tail -f $F | while read n; do echo $((n+1)); sleep 1; done | tee -a $F
