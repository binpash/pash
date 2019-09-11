#!/bin/bash
echo .Enter The Number upto which you want to Print Table: .
read n
i=1
while [ $i -ne 10 ]; do
	i=$(expr $i + 1)
	table=$(expr $i \* $n)
	echo $table
done
