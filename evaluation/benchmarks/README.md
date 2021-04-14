# Experimental Evaluation
Quick Jump: [one-liners](#common-unix-one-liners) | [unix50](#unix-50-from-bell-labs) | [weather analysis](#noaa-weather-analysis) | [web indexing](#wikipedia-web-indexing)

_Most benchmark sets in the evaluation infrastructure include a `input/setup.sh` script for fetching inputs and setting up the experiment appropriately._
See [Running other script]() later.

#### Common Unix one-liners

The one-liner scripts are included in [evaluation/microbenchmarks](../evaluation/microbenchmarks).
The list of scripts (and their correspondence to the names in the paper) are seen below:
 - minimal_grep.sh       # EuroSys: nfa-regex
 - minimal_sort.sh       # EuroSys: sort
 - topn.sh               # EuroSys: top-n
 - wf.sh                 # EuroSys: wf
 - spell.sh              # EuroSys: spell
 - diff.sh               # EuroSys: difference
 - bigrams.sh            # EuroSys: bi-grams
 - set-diff.sh           # EuroSys: set-difference
 - double_sort.sh        # EuroSys: sort-sort
 - shortest_scripts.sh   # EuroSys: shortest-scripts

The inputs that we are going to run them on are defined in
 - `*_env_small.sh` (for the small input)
 - `*_env.sh` (for the large EuroSys eval input, usually 10x larger than the small)

Before running the script we first need to `cd` to the correct directory
  `cd $PASH_TOP/evaluation/eurosys`

The script that runs PaSh on these programs is: [evaluation/eurosys/execute_eurosys_one_liners.sh](../evaluation/eurosys/execute_eurosys_one_liners.sh) 
There are three modes of execution (can be seen by calling the script with the -h flag):

  1. Small inputs | --width 2, 16 | Only full PaSh config
  2. Small inputs | --width 2, 16 | All PaSh configs
  3. Big inputs | -- width 2, 4, 8, 16, 32, 64 | All PaSh configs

The script [evaluation/eurosys/execute_eurosys_one_liners.sh](../evaluation/eurosys/execute_eurosys_one_liners.sh) is based on [evaluation/execute_compile_evaluation_script.sh](../evaluation/execute_compile_evaluation_script.sh) that correctly sets up PaSh for the different configurations.

If you just want to check that PaSh achieves speedups as presented in the paper you can just run 1 with option `-s`.

If you are interested in seeing the improvements by PaSh's runtime primitives (all lines in Figure 9), you can run 2 with option `-m`. 
This should take a couple hours and should validate the trends between different PaSh configurations as shown in Figure 9.

If you want to reproduce the complete results from Figure 9, you need to run 3 with option `-l`.
Note that this should take more than a day to execute.
Also this requires several hundred GBs of free space (due to intermediate inputs, outputs, and buffering).

To plot the results from any of the above experiments, do the following:

```sh
cd $PASH_TOP/compiler
python3 gather_results.py --eurosys2021
```

This will create plots for all invocations of `evaluation/eurosys/execute_eurosys_one_liners.sh`, one for each flag.
The plots are:
* for `-s`: [evaluation/plots/small_tiling_throughput_scaleup.pdf](evaluation/plots/small_tiling_throughput_scaleup.pdf)
* for `-m`: [evaluation/plots/medium_tiling_throughput_scaleup.pdf](../evaluation/plots/medium_tiling_throughput_scaleup.pdf)
* for `-l`: [evaluation/plots/tiling_throughput_scaleup.pdf](../evaluation/plots/tiling_throughput_scaleup.pdf)

Note that `-m` supersedes `-s` but `-l` does not supersede any of the two.

Also note that if you run a script partially, it might end up saving partial results,
therefore having 0 speedups in some points of the plots.

#### Unix50 from Bell Labs

All of the Unix50 pipelines are in [evaluation/unix50/unix50.sh](../evaluation/unix50/unix50.sh).
The inputs of the pipelines are in [evaluation/unix50/](../evaluation/unix50/).

Before running the script we first need to move to the correct directory
  `cd $PASH_TOP/evaluation/eurosys`

The script that runs PaSh on these programs is: [evaluation/eurosys/execute_unix_benchmarks.sh](../evaluation/eurosys/execute_unix_benchmarks.sh) 
There are two modes of execution (can be seen by calling the script with the -h flag):

  1. Small inputs (1GB) | --width 4
  2. Big inputs (10GB) | --width 16 (EuroSys evaluation)

The first one, called with `-s`, uses pash on the unix50 scripts with 1GB input and width 4 
and should be done in less than an hour.
The trend shown in the paper (Fig 10) should be visible in the results from this script.

If you are interested in running the complete evaluation to reproduce Figure 10,
you need to run the script with `-l`. This should take several hours.

To plot the results from any of the above experiments, do the following:

```sh
cd $PASH_TOP/compiler
python3 gather_results.py --eurosys2021
```

This will create plots for both "1GB --width 4" and for "10GB --width 16".

The plots are in:
- for `-s`: [evaluation/plots/unix50_1GB_individual_speedups_4.pdf](../evaluation/plots/unix50_1GB_individual_speedups_4.pdf)
- for `-l`: [evaluation/plots/unix50_10GB_individual_speedups_16.pdf](../evaluation/plots/unix50_10GB_individual_speedups_16.pdf)

Note that the pipelines in the plot are sorted with respect to speedup, and not by their ID.
So the first pipeline does not necessarily correspond to the first pipeline in [evaluation/unix50](../evaluation/unix50).

There are two small differences of these plots compared to Figure 10.
These differences are due to the evolution of PaSh and the refinement of its annotations.
 - First, the first pipeline has higher speedup that 4 and 16 in both cases. This is because
   this pipeline is not very CPU intensive and contains an initial `cat`. PaSh has evolved
   to perform an optimization that removes `cat` occurences that only contain a single file,
   and therefore removes it, improving performance significantly.
 - Second, the slowdown in the last 3 scripts is more significant than the one reported in the paper.
   This is because these scripts contain `tr -d '\n'`, the annotation for which was refined recently due to additional testing.
   The initial annotation for `tr` considered this invocation of `tr` to be stateless while it isn't, 
   since it removes all lines and therefore cannot be parallelized based on lines. The refinement in the annotation
   leads to additional splits to be added after `tr -d '\n'` (since it is non parallelizable pure).
   The issue with these splits is that they do not manage to split the file (since there is only one line)
   leaving the rest of the script to run sequentially.

#### NOAA Weather Analysis

Note that input files that are needed by this script 
are `curl`ed from a server in the local network and therefore
cannot be accessed from elsewhere.

Before running the script we first need to move to the correct directory
  `cd $PASH_TOP/evaluation/eurosys`

The program that we run, described in Section 6.3, can be seen in [evaluation/scripts/max-temp-complete.sh](../evaluation/scripts/max-temp-complete.sh).
It takes as input a sequence of lines each containing a year (e.g. using `seq 2000 2004`).

To run the script with a single year of input use:
  `./execute_max_temp_dish_evaluation.sh -s`

These should take less than 10 minutes.

It runs the script on:
- bash
- pa.sh --width 16

The results are saved in:
- [evaluation/results/max-temp-complete-2000-2000-seq.time](../evaluation/results/max-temp-complete-2000-2000-seq.time)
- [evaluation/results/max-temp-complete-2000-2000-16-pash.time](../evaluation/results/max-temp-complete-2000-2000-16-pash.time)

If you want to run the program with 5 years of input (as is done in Section 6.3)
you need to use the following:
  `./execute_max_temp_dish_evaluation.sh -l`

It should take less than an hour. 
It also runs the script with bash and pash --width 16.

The results are saved in:
- [evaluation/results/max-temp-complete-2000-2004-seq.time](../evaluation/results/max-temp-complete-2000-2004-seq.time)
- [evaluation/results/max-temp-complete-2000-2004-16-pash.time](../evaluation/results/max-temp-complete-2000-2004-16-pash.time)

If you want to separate the preprocessing and processing (as done in Section 6.3)
you need to add the `-e` flag to either 1 or 5 year execution, e.g.:
  `./execute_max_temp_dish_evaluation.sh -l -e`

This runs:
- `evaluation/scripts/max-temp-preprocess.sh`
- `evaluation/scripts/max-temp-process.sh`

with bash, and pash --width 16. It saves results in:
- [evaluation/results/max-temp-preprocess-2000-2000-seq.time](../evaluation/results/max-temp-preprocess-2000-2000-seq.time)
- [evaluation/results/max-temp-preprocess-2000-2000-16-pash.time](../evaluation/results/max-temp-preprocess-2000-2000-16-pash.time)
- [evaluation/results/max-temp-process-2000-2000-seq.time](../evaluation/results/max-temp-process-2000-2000-seq.time)
- [evaluation/results/max-temp-process-2000-2000-16-pash.time](../evaluation/results/max-temp-process-2000-2000-16-pash.time)

and similarly for the large inputs (2000-2004).

Note that PaSh's speedup for the complete script 2000-2004 with width 16
is actually higher than what is reported in the paper since it doesn't
have to write the intermediate files (between preprocessing and processing) to disk.


#### Wikipedia Web Indexing

Note that input files that are needed by this script (complete Wikipedia) 
are saved locally on the server and therefore this program cannot be run from elsewhere.

Before running the script we first need to move to the correct directory
  `cd $PASH_TOP/evaluation/eurosys`

The program that we run, described in Section 6.4, can be seen in [evaluation/scripts/web-index.sh](../evaluation/scripts/web-index.sh).
It requires having set the `$IN`, `$WIKI`, and `$WEB_INDEX_DIR` variables.

To run the script for a 1000 wikipedia links use:
  `./execute_web_index_dish_evaluation.sh -s`

This sets up the required variables and should take less than 5 minutes.
It runs the script with bash, pash --width 2, pash --width 16.

The results are saved in:
- [evaluation/results/web-index-1000-seq.time](../evaluation/results/web-index-1000-seq.time)
- [evaluation/results/web-index-1000-2-pash.time](../evaluation/results/web-index-1000-2-pash.time)
- [evaluation/results/web-index-1000-16-pash.time](../evaluation/results/web-index-1000-16-pash.time)

If you want to run with the EuroSys evaluation inputs (100k links), use:
  `./execute_web_index_dish_evaluation.sh -l`

This should take a couple hours and the results are saved in:
- [evaluation/results/web-index-100000-seq.time](../evaluation/results/web-index-100000-seq.time)
- [evaluation/results/web-index-100000-2-pash.time](../evaluation/results/web-index-100000-2-pash.time)
- [evaluation/results/web-index-100000-16-pash.time](../evaluation/results/web-index-100000-16-pash.time)
