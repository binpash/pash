#!/bin/bash
uniq -c | sort -nr | head -n 1 | fmt -w1 | sed 1d
