#!/bin/bash

## Setup PASH_TOP if it is not set
export PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}

## Setup inputs if they are not there already
echo "Downloading/setting up input..."
cd "$PASH_TOP/evaluation/benchmarks/oneliners/input"
# remove old input
./setup.sh -c
./setup.sh --full

## Set parallelism level
## Warning: these two have to mirror each other
export JOBS=16
export IN=$(echo in{00..15})
export IN=$(echo $IN)

## Setup input and dictionary
export IN_FILE="$PASH_TOP/evaluation/benchmarks/oneliners/input/1G.txt"
export dict="$PASH_TOP/evaluation/benchmarks/oneliners/input/dict.txt"

cd "$PASH_TOP/evaluation/micro/gnu-parallel"
## Split the files
split -n l/$JOBS -d "$IN_FILE" in

echo "Running spell with bash..."
time ./spell.sh > spell.bash.out

run_spell_gnu_parallel()
{
    echo "Running spell with gnu-parallel with block_size=$BLOCK_SIZE..."
    time ./spell.par.sh > spell.parallel.$BLOCK_SIZE.out
    echo "Comparing bash and gnu-parallel output"
    diff spell.bash.out spell.parallel.$BLOCK_SIZE.out | head
}

export BLOCK_SIZE=50K
run_spell_gnu_parallel

export BLOCK_SIZE=250K
run_spell_gnu_parallel

export BLOCK_SIZE=250M
run_spell_gnu_parallel

echo "Running spell with pash..."
time "$PASH_TOP/pa.sh" --r_split_batch_size 1000000 --r_split --dgsh_tee -w $JOBS ./spell.sh > spell.pash.out
echo "Comparing bash and pash output"
diff spell.bash.out spell.pash.out | head



## Setup input for set-diff
export IN_FILE="$PASH_TOP/evaluation/benchmarks/oneliners/input/3G.txt"
split -n l/$JOBS -d "$IN_FILE" in

echo "Running set-diff with bash..."
time ./set-diff.sh > set-diff.bash.out

echo "Running set-diff with gnu-parallel with block_size=$BLOCK_SIZE..."
time ./set-diff.par.sh > set-diff.parallel.$BLOCK_SIZE.out
echo "Comparing bash and gnu-parallel output"
diff set-diff.bash.out set-diff.parallel.$BLOCK_SIZE.out | head

echo "Running set-diff with pash..."
time "$PASH_TOP/pa.sh" --r_split_batch_size 1000000 --r_split --dgsh_tee -w $JOBS ./set-diff.sh > set-diff.pash.out
echo "Comparing bash and pash output"
diff set-diff.bash.out set-diff.pash.out | head
