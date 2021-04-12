#### Section 6.5: Further Micro-benchmarks

To run the comparison with sort --parallel, just use [evaluation/eurosys/execute_baseline_sort.sh](../evaluation/eurosys/execute_baseline_sort.sh)

Before running the script we first need to move to the correct directory
  `cd $PASH_TOP/evaluation/eurosys`

There are two modes of execution:
1. option: -s Small input | --width 2, 16
2. option: -l Big input | -- width 2, 4, 8, 16, 32, 64

Note that this script executes sort --parallel with double the value of --width
since we noticed that it grows slightly slower (as shown in the Figure in Section 6.5).

_This script throws a warning that is expected: `Env file: .../evaluation/microbenchmarks/sort_env_small.sh could not be read.` The warning is expected and can be safely ignored._
