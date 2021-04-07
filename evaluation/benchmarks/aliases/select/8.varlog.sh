#!/bin/bash 
set -e
# 1308; or line above, w/ -vE
# Doesn't do much
find /var/log -type f -exec file {} \; | grep 'text' | cut -d' ' -f1 | sed -e's/:$//g' | grep -v '[0-9]$' | xargs tail 
