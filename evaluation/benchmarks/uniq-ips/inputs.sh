#!/bin/bash
cd $(dirname "$0")

N=200000000 # 200M

mkdir -p inputs
./gen_data.py "$N" > inputs/logs-popcount-org-large.txt

"$PASH_TOP/scripts/append_nl_if_not.sh" inputs/logs-popcount-org-large.txt