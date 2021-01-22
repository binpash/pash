#!/bin/bash

## Things to note:
## - gen_big_files has to run together with install
## - PaSh is storage hungry!!! So we need to give reviewers access to machines with tons of storage
##   + Either give them access to 3 EC2 instances with 64 instances and 1TB storage (they might cost a lot)
##   + Or give all of them access to deathstar
##   + Both of the above should provide a `screen` and instructions about how to use it.

## TODO: Copy the plots (Unix50 1GB, 10GB and the one-liners in a different dir)



## TODO: By default run all small
##       In the video, we need to explain all options, and what do you get for all options.

## TODO: Set up $PASH_TOP in the beginning or run the install script.

## TODO: Note that the execution times reported in the paper are in the log (and not the general time commands).

echo "This script runs the whole EuroSys 2021 PaSh evaluation"

## TODO: Get the provided execution level
## TODO: Can also have an only-plot option  


echo ""
echo "Section 6.1: Common Unix One-liners"
## TODO: Run the one-liners evaluation

## TODO: Run the plotting 
## TODO: Also save aggregates (avg, etc) in a file
## TODO: Depending on flag, echo where the results and plot are

## Note that input files that are used as inputs for this script are generated 
## using the `gen*` scripts in `evaluation/scripts/input/`.
## ```
## cd $PASH_TOP/evaluation/scripts/input/
## ./gen.sh
## ./gen.sh # Warning: This requires more than 100GB of space.
## ```
##
## If you just want to run the scripts with small inputs (the main conclusions still hold)
## you only need to run `./gen.sh`.
##
## The one-liner scripts are included in `evaluation/microbenchmarks`
## The list of scripts (and their correspondence to the names in the paper) are seen below:
##  - minimal_grep.sh       # EuroSys: nfa-regex
##  - minimal_sort.sh       # EuroSys: sort
##  - topn.sh               # EuroSys: top-n
##  - wf.sh                 # EuroSys: wf
##  - spell.sh              # EuroSys: spell
##  - diff.sh               # EuroSys: difference
##  - bigrams.sh            # EuroSys: bi-grams
##  - set-diff.sh           # EuroSys: set-difference
##  - double_sort.sh        # EuroSys: sort-sort
##  - shortest_scripts.sh   # EuroSys: shortest-scripts
##
## The inputs that we are going to run them on are defined in
##  - *_env_small.sh (for the small input)
##  - *_env.sh (for the large EuroSys eval input, usually 10x larger than the small)
##
## The script that runs PaSh on these programs is: `evaluation/eurosys/execute_eurosys_one_liners.sh` 
## There are three modes of execution (can be seen by calling the script with the -h flag):
##   1. Small inputs | --width 2, 16 | Only full PaSh config
##   2. Small inputs | --width 2, 16 | All PaSh configs
##   3. Big inputs | -- width 2, 4, 8, 16, 32, 64 | All PaSh configs
##
## If you just want to check that PaSh achieves speedups as presented in the paper
## you can just run 1 with option `-s`.
##
## If you are interested in seeing the improvements by PaSh's runtime primitives
## (all lines in Figure 9), you can run 2 with option `-m`. 
## This should take a couple hours and should validate the trends between different PaSh
## configurations as shown in Figure 9.
##
## If you want to reproduce the complete results from Figure 9, you need to run 3 with option `-l`.
## Note that this should take more than a day to execute.
## Also this requires several hundred GBs of free space (due to intermediate inputs, outputs, and buffering).
##
## To plot the results from any of the above experiments, do the following:
## ```
## cd $PASH_TOP/compiler
## python3 gather_results.py
## ```
##
## This will create plots for all invocations of 
## `evaluation/eurosys/execute_eurosys_one_liners.sh`, one for each flag.
##
## The plots are:
## - for `-s`: evaluation/plots/small_tiling_throughput_scaleup.pdf
## - for `-m`: evaluation/plots/medium_tiling_throughput_scaleup.pdf
## - for `-l`: evaluation/plots/tiling_throughput_scaleup.pdf
##
## Note that `-m` supersedes `-s` but `-l` does not supersede any of the two.
##
## Also note that if you run a script partially, it might end up saving partial results,
## therefore having 0 speedups in some points of the plots.

echo ""
echo "Section 6.2: Unix50 from Bell Labs"
## TODO: Run the unix50

## TODO: Run the plotting 
## TODO: Also save aggregates (avg, etc) in a file
## TODO: Depending on flag, echo where the results and plot are

## TODO: Inputs of max-temp have to be accessible...

echo ""
echo "Section 6.3: Use Case: NOAA Weather Analysis"

## Note that input files that are needed by this script 
## are `curl`ed from a server in the local network and therefore
## cannot be accessed from elsewhere.
##
## Before running the script we first need to move to the correct directory
##   `cd $PASH_TOP/evaluation/eurosys`
##
## The program that we run, described in Section 6.3, can be seen in `evaluation/scripts/max-temp-complete.sh`.
## It takes as input a sequence of lines each containing a year (e.g. using `seq 2000 2004`).
##
## To run the script with a single year of input use:
##   `./execute_max_temp_dish_evaluation.sh -s`
##
## These should take less than 10 minutes.
##
## It runs the script on:
## - bash
## - pa.sh --width 16
##
## The results are saved in:
## - `evaluation/results/max-temp-complete-2000-2000-seq.time`
## - `evaluation/results/max-temp-complete-2000-2000-16-pash.time`
##
## If you want to run the program with 5 years of input (as is done in Section 6.3)
## you need to use the following:
##   `./execute_max_temp_dish_evaluation.sh -l`
##
## It should take less than an hour. 
## It also runs the script with bash and pash --width 16.
##
## The results are saved in:
## - `evaluation/results/max-temp-complete-2000-2004-seq.time`
## - `evaluation/results/max-temp-complete-2000-2004-16-pash.time`
##
## If you want to separate the preprocessing and processing (as done in Section 6.3)
## you need to add the `-e` flag to either 1 or 5 year execution, e.g.:
##   `./execute_max_temp_dish_evaluation.sh -l -e`
##
## This runs:
## - `evaluation/scripts/max-temp-preprocess.sh`
## - `evaluation/scripts/max-temp-process.sh`
##
## with bash, and pash --width 16. It saves results in:
## - `evaluation/results/max-temp-preprocess-2000-2000-seq.time`
## - `evaluation/results/max-temp-preprocess-2000-2000-16-pash.time`
## - `evaluation/results/max-temp-process-2000-2000-seq.time`
## - `evaluation/results/max-temp-process-2000-2000-16-pash.time`
##
## and similarly for the large inputs (2000-2004).
##
## Note that PaSh's speedup for the complete script 2000-2004 with width 16
## is actually higher than what is reported in the paper since it doesn't
## have to write the intermediate files (between preprocessing and processing) to disk.
##

echo ""
echo "Section 6.4: Use Case: Wikipedia Web Indexing"

## Note that input files that are needed by this script (complete Wikipedia) 
## are saved locally on the server and therefore this program cannot be run from elsewhere.
##
## Before running the script we first need to move to the correct directory
##   `cd $PASH_TOP/evaluation/eurosys`
##
## The program that we run, described in Section 6.4, can be seen in `evaluation/scripts/web-index.sh`.
## It requires having set the `$IN`, `$WIKI`, and `$WEB_INDEX_DIR` variables.
##
## To run the script for a 1000 wikipedia links use:
##   `./execute_web_index_dish_evaluation.sh -s`
##
## This sets up the required variables and should take less than 5 minutes.
## It runs the script with bash, pash --width 2, pash --width 16.
##
## The results are saved in:
## - `evaluation/results/web-index-1000-seq.time`
## - `evaluation/results/web-index-1000-2-pash.time`
## - `evaluation/results/web-index-1000-16-pash.time`
##
## If you want to run with the EuroSys evaluation inputs (100k links), use:
##   `./execute_web_index_dish_evaluation.sh -l`
##
## This should take a couple hours and the results are saved in:
## - `evaluation/results/web-index-100000-seq.time`
## - `evaluation/results/web-index-100000-2-pash.time`
## - `evaluation/results/web-index-100000-16-pash.time`

echo ""
echo "Section 6.5: Further Micro-benchmarks"

## TODO: Run the parallel sort

## TODO(@nikos): Run and explain the GNU Parallel
