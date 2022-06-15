#!/bin/bash
uniq -c | sort -nr | head -n 1 | tr -s ' ' '\n' | tail -n 1
