#!/bin/bash

## Note: Needs to be run on a big git repository to make sense (maybe linux)

## Initialize the necessary temporary files
file1=$(mktemp)

forder()
{
	sort |
	uniq -c |
	sort -rn
}


git log --format="%an:%ad" --date=default "$@" >"$file1"

echo "Authors ordered by number of commits"
# Order by frequency
awk -F: '{print $1}' <"$file1" | forder

echo "Days ordered by number of commits"
# Order by frequency
awk -F: '{print substr($2, 1, 3)}' <"$file1" | forder
