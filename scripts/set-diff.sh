#!/bin/bash
# Show diff between different groups of commands
f1=${1:-posix}
f2=${2:-coreutils}
p="../stats/c_stats"
comm -3 <(cut -d ' ' -f 1 $p/$f1.txt | sort ) <( cut -d ' ' -f 1 $p/$f2.txt | sort)
