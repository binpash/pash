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

## TODO: Ensure that all input files (e.g. x100 for spell) are present.

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
## TODO: Run the weather analysis max-temp

## TODO: Also save results in files and report


## TODO: Inputs of web-indexing have to be accessible...


echo ""
echo "Section 6.4: Use Case: Wikipedia Web Indexing"
## TODO: Run the web-index script

## Also save the results in files and report
