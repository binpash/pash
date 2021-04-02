#!/bin/bash
# First command is almost always a generator
set -e
IN=${PASH_TOP}/evaluation/benchmarks/aliases/meta/
#FIXME define a complex expression
COMPLEX=""
# linux git
cd ${IN}/linux
git ls-tree --name-only -z -r HEAD | grep -z -Z -E '\.(cc|h|cpp|hpp|c|txt|java)$' | xargs -0 -n1 git blame --line-porcelain | grep ${COMPLEX}

