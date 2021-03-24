#!/bin/bash
set -e

# time: print real in seconds, to simplify parsing
TIMEFORMAT="%3R" # %3U %3S"

if [[ -z "$PASH_TOP" ]]; then
  echo "Need to provide PASH_TOP, possibly $(git rev-parse --show-toplevel)" 1>&2
  exit 1
fi

echo executing one-liners
export IN=${IN:-$PASH_TOP/evaluation/benchmarks/oneliners/input/1M.txt}
cd oneliners/
echo nfa-regex.sh:        $({ time ./nfa-regex.sh > /dev/null; } 2>&1)
echo sort.sh:             $({ time ./sort.sh > /dev/null; } 2>&1)
echo top-n.sh:            $({ time ./top-n.sh > /dev/null; } 2>&1)
echo wf.sh:               $({ time ./wf.sh > /dev/null; } 2>&1)
echo spell.sh:            $({ time ./spell.sh > /dev/null; } 2>&1)
echo bi-grams.sh:         $({ time ./bi-grams.sh > /dev/null; } 2>&1)
echo diff.sh:             $({ time ./diff.sh > /dev/null; } 2>&1)
echo set-diff.sh:         $({ time ./set-diff.sh > /dev/null; } 2>&1)
echo shortest-scripts.sh: $({ time ./shortest-scripts.sh > /dev/null; } 2>&1)
echo sort-sort.sh:        $({ time ./sort-sort.sh > /dev/null; } 2>&1)

