#!/bin/bash
# Clever trick that uses the /dev/fd/xx pseudo-file system
# https://stackoverflow.com/questions/40244/how-to-make-a-pipe-loop-in-bash

# MMG 2022-06-30 the `function` kw is a bash-ism; leaving it in to not disrupt what gets optimized in previous evaluations
function calc() {
  # calculate sum of squares of numbers 0,..,10

  sum=0
  for ((i=0; i<10; i++)); do
    echo $i                   # "request" the square of i

    read ii                   # read the square of i
    echo "got $ii" >&2          # debug message

    let sum=$sum+$ii
  done

  echo "sum $sum" >&2           # output result to stderr
}

function square() {
  # square numbers

  read j                         # receive first "request"
  while [ "$j" != "" ]; do
    let jj=$j*$j
    echo "square($j) = $jj" >&2  # debug message

    echo $jj                     # send square

    read j                       # receive next "request"
  done
}

read | { calc | square; } >/dev/fd/0
