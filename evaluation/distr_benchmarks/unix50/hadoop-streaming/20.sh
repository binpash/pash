#!/bin/bash
cat $1 | grep '(' | cut -d '(' -f 2 | cut -d ')' -f 1
