#!/bin/bash
cat $1 | cut -c 1-2 | tr -d '\n'
