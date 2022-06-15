#!/bin/bash
cat $1 | grep 'UNIX' | cut -f 1
