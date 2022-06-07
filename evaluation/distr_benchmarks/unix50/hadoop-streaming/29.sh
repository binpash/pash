#!/bin/bash
cat $1 | sed 2d | sed 2d | tr -c '[A-Z]' '\n' | tr -d '\n'
