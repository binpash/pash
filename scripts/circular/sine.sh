#!/bin/bash
[ -f /etc/resolv.conf ] || ( echo 1 >temp.txt )
tail -f temp.txt | while read n; do echo "1+s(3*$n)" | bc -l; sleep 1; done | tee -a temp.txt
