#!/bin/bash

echo 'GNU Coreutils ('$(cat coreutils-summary.txt | wc -l | awk '{$1=$1};1') 'commands):'
echo '  S:' $(cat coreutils-summary.txt | grep ' S ' | wc -l)
echo '  P:' $(cat coreutils-summary.txt | grep ' P ' | wc -l)
echo '  N:' $(cat coreutils-summary.txt | grep ' N ' | wc -l)
echo '  E:' $(cat coreutils-summary.txt | grep ' E ' | wc -l)

echo 'POSIX ('$(( $(cat posix-summary.txt | wc -l | awk '{$1=$1};1') + $(cat ../c_stats/posix.txt | grep -v Mandatory | wc -l) )) 'commands):'
echo '  S:' $(cat posix-summary.txt | grep ' S ' | wc -l)
echo '  P:' $(cat posix-summary.txt | grep ' P ' | wc -l)
echo '  N:' $(cat posix-summary.txt | grep ' N ' | wc -l)
echo '  E:' $(( $(cat posix-summary.txt | grep ' E ' | wc -l) + $(cat ../c_stats/posix.txt | grep -v Mandatory | wc -l) ))

