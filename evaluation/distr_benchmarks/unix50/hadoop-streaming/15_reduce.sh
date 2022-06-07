#!/bin/bash
i=0
while read line
do
  ((i=i+$line))
done
echo $i
