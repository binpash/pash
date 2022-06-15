#!/bin/bash
cat $1 | grep 'Bell' | cut -f 2
