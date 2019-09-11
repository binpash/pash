#!/bin/bash
x=0
y=1
i=2
while true ; do
	i=$(expr $i + 1)
	z=$(expr $x + $y)
	echo -n "$z "
	x=$y
	y=$z
	sleep .5
done
