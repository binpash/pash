#!/bin/bash
cat $1 | cut -c 1-1 | tr -d '\n'
