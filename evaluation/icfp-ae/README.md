
Before going further, please skim the [online tutorial](https://github.com/andromeda/pash/blob/icfp-ae/docs/tutorial.md) to get a sense of where to look for details if something goes wrong. For reviewers setting up PaSh on own infrastructure, the part of the tutorial outlining [installation instructions](https://github.com/andromeda/pash/blob/icfp-ae/docs/tutorial.md#installation) should be particularly relevant.

# Claims

The paper outlines four main contributions, of which only the last is related to the artifact (p.2 of the submitted PDF):

> * Results: It reimplements the optimizing compiler of PaSh and presents experimental results that evaluate the speedup afforded by a specific dataflow transformation (§7).

We claim all three badges for this contribution: Available, Functional, and Reproducible.

##  "Artifact Available" badge: 

Our artifact is available on GitHub, including all the changes described in the paper, permanently and publicly. The re-implementation discussed in the paper, including the bug fixes, includes several commits---the majority of which are included in [Pull Request #83](https://github.com/andromeda/pash/commit/94b09e71316e8a0b10e0b6450a0b2953a04a71df). This PR and other fixes have been merged into PaSh, and are permanantely and publicly available in PaSh's GitHub repository.

##  "Artifact Functional" badge: 

PaSh, including the changes described in the paper, is fully functional. It includes all key components described in the paper, automated scripts, and input data for to run the experiments described in the paper. It also includes extensive documentation: [Introduction](https://github.com/andromeda/pash/blob/icfp-ae/docs/tutorial.md#introduction) | [Installation](https://github.com/andromeda/pash/blob/icfp-ae/docs/tutorial.md#installation) | [Running Scripts](https://github.com/andromeda/pash/blob/icfp-ae/docs/tutorial.md#running-scripts) | [What Next?](https://github.com/andromeda/pash/blob/icfp-ae/docs/tutorial.md#what-next)

Each of the main sets of evaluation benchmarks is in its own directory. For benchmarks in §7.1:

* The GNU `parallel` case study can be found in [ ](https://github.com/andromeda/pash/blob/icfp-ae/evaluation/icfp-ae/gnu-parallel-case-study.sh)

For benchmarks in §7.2:

* Expert Pipelines can be found in [./evaluation/benchmarks/oneliners](https://github.com/andromeda/pash/tree/icfp-ae/evaluation/benchmarks/oneliners)
* Unix50 Pipelines can be found in [./evaluation/benchmarks/unix50](https://github.com/andromeda/pash/tree/icfp-ae/evaluation/benchmarks/unix50)
* COVID-19 Mass-Transit Analysis Pipelines can be found in [./evaluation/benchmarks/unix50](https://github.com/andromeda/pash/tree/icfp-ae/evaluation/benchmarks/analytics-mts)

Note that each one of the last three benchmark sets has its own `./input/setup.sh` in the directory of the benchmark for setting up its inputs.

##  "Artifact Reproduced" badge: 

To reproduce the results presented in the paper:
1. Log into the server, as outlined in the secret gist, _OR_ Follow the [installation instructions](https://github.com/andromeda/pash/blob/icfp-ae/docs/tutorial.md#installation) to set up PaSh on your infrastructure
2. Make sure you have the `PASH_TOP` environment variable set up pointing to the top-level directory in the repo: `echo $PASH_TOP`.
3. Run the following commands:

```sh
$PASH_TOP/evaluation/icfp-ae/gnu-parallel-case-study.sh
$PASH_TOP/evaluation/icfp-ae/oneliners.sh
$PASH_TOP/evaluation/icfp-ae/unix50.sh
$PASH_TOP/evaluation/icfp-ae/analytics-mts.sh
```

Passing no flags to these scripts will execute smaller inputs (in under one hour), so that the reviewers can quickly confirm that everything works as expected.
Passing the `--full` flag on all benchmarks will run with the full input set, as reported in the paper.
(Note: these inputs have already been set up on `deathstar`, the machine used for , but will take additional time to download and set up on other environments and will require about 300GB).

The results include the sequential baseline (running `bash`), `pash` without the `cat`-`split` transformation, and full `pash` performing all transformations. The scripts are configured to apply `pash` with a parallelism `--width` of `16`.

