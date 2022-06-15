#!/bin/bash
cat $1 | sed 1d | grep 'Bell' | cut -f 2
