#!/bin/bash
cat $IN7 | cut -f 1 | grep 'AT&T' | wc -l
