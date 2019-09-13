#!/bin/bash
# https://unix.stackexchange.com/questions/193441/how-can-i-implement-a-circular-flow-of-data-among-interconnected-commands
echo 1 >temp.txt
tail -f temp.txt | while read n; do echo $((n+1)); sleep 1; done | tee -a temp.txt
