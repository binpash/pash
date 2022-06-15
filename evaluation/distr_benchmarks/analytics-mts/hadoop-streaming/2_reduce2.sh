#!/bin/bash
uniq -c | sort -k1n | awk -v OFS="\t" "{print \$2,\$1}"
