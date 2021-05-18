#!/bin/bash

## Setup PASH_TOP if it is not set
export PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}

## Setup inputs if they are not there already
cd "$PASH_TOP/evaluation/benchmarks/oneliners/input"
./setup.sh

## Set parallelism level
## Warning: these two have to mirror each other
export JOBS=16
export IN=$(echo in{00..15})
export IN=$(echo $IN)

## Setup input and dictionary
export IN_FILE="$PASH_TOP/evaluation/benchmarks/oneliners/input/10M.txt"
export dict="$PASH_TOP/evaluation/benchmarks/oneliners/input/dict.txt"

cd "$PASH_TOP/evaluation/micro/gnu-parallel"
## Split the files
split -n l/$JOBS -d "$IN_FILE" in

echo "Running bash..."
time ./spell.sh > spell.bash.out

run_spell_gnu_parallel()
{
    echo "Running gnu-parallel with block_size=$BLOCK_SIZE..."
    time ./spell.par.sh > spell.parallel.$BLOCK_SIZE.out
    echo "Comparing bash and gnu-parallel output"
    diff spell.bash.out spell.parallel.$BLOCK_SIZE.out    
}

export BLOCK_SIZE=50K
run_spell_gnu_parallel

export BLOCK_SIZE=250K
run_spell_gnu_parallel

export BLOCK_SIZE=250M
run_spell_gnu_parallel

echo "Running pash..."
time "$PASH_TOP/pa.sh" -w $JOBS ./spell.sh > spell.pash.out
echo "Comparing bash and pash output"
diff spell.bash.out spell.pash.out
