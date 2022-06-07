#!/bin/bash
cat $1 | cut -f 2 | cut -d ' ' -f 1
