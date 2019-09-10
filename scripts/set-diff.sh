#!/bin/bash
p="../linux_command_statistics/c_stats"
echo $1 $2 $3
comm  -3 <(cut -d ' ' -f 1 $p/$1.txt | sort ) <( cut -d ' ' -f 1 $p/$2.txt | sort)
