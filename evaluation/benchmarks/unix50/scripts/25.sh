#!/bin/bash

# 9.3: animal that used to decorate the Unix room
cat $IN | cut -c 1-2 | tr -d '\n'
