#!/bin/bash

# To build and sign an SGX enclave, this script extracts the executable's symbol
# table and calculates its SHA256 hashsum.

# Require
# Data: /usr/lib/libz3.so

IN=/usr/lib/libz3.so
OUT=./output/out.txt

readelf -x .symtab $IN |
  tail -n +3 |
  head -n -1 |                # next three implement `awk '{print $2$3$4$5}'`
  sed 's/^[[:space:]]*//' | 
  cut -d ' ' -f2-5 |
  tr -d ' ' |
  tr -d "\n" |
  xxd -r -p |
  sha256sum > $OUT
