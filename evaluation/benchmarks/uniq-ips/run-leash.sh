#!/bin/bash

cd "$(dirname "$0")" || exit 1

WIDTH=32
time IN=uniq-ips/inputs/logs-popcount-org.txt OUT=uniq-ips/outputs/uniq-ips.sh:$WIDTH:hybrid: "$PASH_TOP"/pa.sh -w "$WIDTH" scripts/uniq-ips.sh --serverless_exec