#!/bin/bash

oldifs=$IFS
IFS=$(echo -e "\t")
for f in `ls *`; do
    echo $f
done
IFS=$oldifs
