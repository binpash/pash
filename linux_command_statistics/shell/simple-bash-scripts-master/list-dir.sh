#!/bin/bash

list=($(ls))

for f in "${list[@]}";do
    echo $f
done