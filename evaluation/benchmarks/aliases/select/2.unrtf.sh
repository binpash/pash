#!/bin/bash

#tag: rtf-to-txt
set -e

IN=$PASH_TOP/evaluation/benchmarks/aliases/meta/rtf
OUT=$PASH_TOP/evaluation/benchmarks/aliases/meta/out
find $IN -name '*.rtf' | xargs -I {} unrtf {} --text > /dev/null
