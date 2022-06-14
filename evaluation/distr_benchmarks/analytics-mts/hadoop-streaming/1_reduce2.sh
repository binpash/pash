#!/bin/bash
uniq -c | awk -v OFS="\t" "{print \$2,\$1}"
