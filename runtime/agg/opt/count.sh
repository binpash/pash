#!/usr/bin/env bash
awk '{ count[$2] += $1 } END { for(e in count) print count[e], e }' "$@"
