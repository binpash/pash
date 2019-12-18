#!/bin/bash
name=$1
path=$2
tar -czvf $name.tar.gz $path 
gpg -c $name.tar.gz
rm -rf $name.tar.gz
