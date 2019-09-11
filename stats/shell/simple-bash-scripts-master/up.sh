#!/bin/bash

LEVEL=$1
for ((i = 0; i < LEVEL; i++)); do
	echo $i
	CDIR=../$CDIR
done
cd $CDIR
echo "You are in: "$PWD
sh=$(which $SHELL)
exec $sh
