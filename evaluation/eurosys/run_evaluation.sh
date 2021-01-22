#!/bin/bash

## BIG Issues:
## - Inputs for max-temp have to be accessible
## - Input for web-indexing has to be copied
## - gen_big_files has to run together with install
## - PaSh is storage hungry!!! So we need to give reviewers access to machines with tons of storage
##   + Either give them access to 3 EC2 instances with 64 instances and 1TB storage (they might cost a lot)
##   + Or give all of them access to deathstar
##   + Both of the above should provide a `screen` and instructions about how to use it.

## TODO: Copy the plots (Unix50 1GB, 10GB and the one-liners in a different dir)



## TODO: By default run all small
##       In the video, we need to explain all options, and what do you get for all options.



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
## are curled from a server in the local network and therefore
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
## TODO: Run the web-index script

## Also save the results in files and report
