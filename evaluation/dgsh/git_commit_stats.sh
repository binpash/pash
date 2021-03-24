#!/bin/bash
# Adapted from the DGSH
# https://www.spinellis.gr/sw/dgsh/#commit-stats

# Setup:
# git clone https://github.com/torvalds/linux && cd linux

# A function the original script uses, in case we want to use it
forder()
{
  sort |
  uniq -c |
  sort -rn
}

mkfifo a b
cat a | awk -F: '{print $1}' | sort | uniq -c | sort -rn | { echo Authors ordered by number of commits ; cat ; } &
cat b | awk -F: '{print substr($2, 1, 3)}' | sort | uniq -c | sort -rn | { echo Days ordered by number of commits ; cat ; } &
git log --format="%an:%ad" --date=default "$@" | tee a b >/dev/null &
rm a b
