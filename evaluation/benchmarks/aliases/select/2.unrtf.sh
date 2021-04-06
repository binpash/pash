#!/bin/bash

#tag: rtf-to-txt
set -e
IN=${RTF:-$PASH_TOP/evaluation/benchmarks/aliases/meta/rtf}
OUT=${OUT:-PASH_TOP/evaluation/benchmarks/aliases/meta/out}
find $IN -name '*.rtf' | xargs -I {} unrtf {} --text > /dev/null
