#!/usr/bin/env bash

## Multiplies its stdin by -m times

multiply_factor=1

## TODO: Implement
size_limit=0

while getopts 'm:l' opt; do
    case $opt in
        m) multiply_factor=$OPTARG ;;
        l) echo "Option -l not implemented yet"
           exit 1 ;;
        *) echo 'Error in command line parsing' >&2
           exit 1
    esac
done
shift "$(( OPTIND - 1 ))"

# echo "Mult by: $multiply_factor"

temp_file="$(mktemp -u)"

cat > temp_file

for (( i = 0; i < $multiply_factor; i++ )); do
  cat temp_file
done

rm temp_file
