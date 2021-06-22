
Before going further, please skim the [online tutorial](https://github.com/andromeda/pash/blob/icfp-ae/docs/tutorial.md) to get a sense of where to look for details if something goes wrong. For reviewers setting up PaSh on own infrastructure, the part of the tutorial outlining [installation instructions](https://github.com/andromeda/pash/blob/icfp-ae/docs/tutorial.md#installation) should be particularly relevant.

# Claims

The paper outlines four main contributions, of which only the last is related to the artifact (p.2 of the submitted PDF):

> * Results: It reimplements the optimizing compiler of PaSh and presents experimental results that evaluate the speedup afforded by a specific dataflow transformation (§7).

We claim all three badges for this contribution: Available, Functional, and Reproducible.

##  "Artifact Available" badge: 

Our artifact is available on GitHub, including all the changes described in the paper, permanently and publicly. The re-implementation discussed in the paper, including the bug fixes, includes several commits---the majority of which are included in [Pull Request #83](https://github.com/andromeda/pash/commit/94b09e71316e8a0b10e0b6450a0b2953a04a71df). This PR and other fixes have been merged into PaSh, and are permanantely and publicly available in PaSh's GitHub repository.

For archival purposes, we have also uploaded a [quemu VM image on Zenodo](https://zenodo.org/record/4776838). We do not recommend using this VM for reproducing the exact paper results though, due to performance overheads associated with the VM.

##  "Artifact Functional" badge: 

PaSh, including the changes described in the paper, is fully functional. It includes all key components described in the paper, automated scripts, and input data for to run the experiments described in the paper. It also includes extensive documentation: [Introduction](https://github.com/andromeda/pash/blob/icfp-ae/docs/tutorial.md#introduction) | [Installation](https://github.com/andromeda/pash/blob/icfp-ae/docs/tutorial.md#installation) | [Running Scripts](https://github.com/andromeda/pash/blob/icfp-ae/docs/tutorial.md#running-scripts) | [What Next?](https://github.com/andromeda/pash/blob/icfp-ae/docs/tutorial.md#what-next)

Each of the main sets of evaluation benchmarks is in its own directory. For benchmarks in §7.1:

* The GNU `parallel` case study can be found in [./evaluation/icfp-ae/gnu-parallel-case-study.sh](https://github.com/andromeda/pash/blob/icfp-ae/evaluation/icfp-ae/gnu-parallel-case-study.sh)

For benchmarks in §7.2:

* Expert Pipelines can be found in [./evaluation/benchmarks/oneliners](https://github.com/andromeda/pash/tree/icfp-ae/evaluation/benchmarks/oneliners)
* Unix50 Pipelines can be found in [./evaluation/benchmarks/unix50](https://github.com/andromeda/pash/tree/icfp-ae/evaluation/benchmarks/unix50)
* COVID-19 Mass-Transit Analysis Pipelines can be found in [./evaluation/benchmarks/analytics-mts](https://github.com/andromeda/pash/tree/icfp-ae/evaluation/benchmarks/analytics-mts)

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

The results include the sequential baseline (running `bash`), `pash` without the `cat`-`split` transformation, and full `pash` performing all transformations. The scripts are configured to apply `pash` with a parallelism `--width` of `16`.

### Notes 

- The input sizes for `oneliners` have been scaled down (x3 or x10) with respect to the ones shown in the paper so that they execute in a reasonable amount of time.

- The input size and block size for the `spell` part of the gnu parallel evaluation is also slightly modified so that it executes in a reasonable amount of time (input size was reduced from 3G to 1G and block size 10K was modified to 50K).

- The unix50 scripts (`1.sh` - `36.sh`) contain two scripts that are empty (`22.sh` and `27.sh`) that are not shown in the paper plot.

- Note on hardware and software requirements: Running these scripts will require significant disk space (>100GB), it will take considerable time (several hours), and will need many CPUs (ideally more than 32). AEC reviewers are provided access to `deathstar`, the machine used to gather the results reported in the paper, which already meets these hardware and software requirements.

- Our claims and results are not particular to specific command versions, but [pkg_versions.txt](../scripts/pkg_versions.txt) contains versions of requirements on our experimental infrastructure for future reproducibility.


### Expected Results

We collect the output of a run of all experiments for future reference:

#### GNU Parallel

```
Downloading/setting up input...
Generating full-size inputs
Running spell with bash...

real    8m30.196s
user    9m3.771s
sys     2m34.787s
Running spell with gnu-parallel with block_size=50K...

real    2m48.019s
user    11m51.651s
sys     6m28.305s
Comparing bash and gnu-parallel output
Running spell with gnu-parallel with block_size=250K...

real    1m10.121s
user    9m37.386s
sys     3m11.701s
Comparing bash and gnu-parallel output
Running spell with gnu-parallel with block_size=250M...

real    0m52.812s
user    9m0.701s
sys     2m24.801s
Comparing bash and gnu-parallel output
Running spell with pash...

real    1m2.798s
user    11m3.757s
sys     3m38.919s
Comparing bash and pash output
Running set-diff with bash...

real    17m6.116s
user    18m19.996s
sys     0m25.974s
Running set-diff with gnu-parallel with block_size=250M...

real    2m40.809s
user    18m30.232s
sys     0m38.266s
Comparing bash and gnu-parallel output
Running set-diff with pash...

real    2m20.964s
user    22m54.571s
sys     2m33.453s
Comparing bash and pash output
```

#### Oneliners

```
nfa-regex.sh:                  478.333  nfa-regex.sh:                  31.031   nfa-regex.sh:                  30.442
sort.sh:                       376.471  sort.sh:                       75.852   sort.sh:                       75.187
top-n.sh:                      412.159  top-n.sh:                      101.036  top-n.sh:                      83.411
wf.sh:                         402.054  wf.sh:                         109.512  wf.sh:                         88.821
spell.sh:                      479.021  spell.sh:                      152.901  spell.sh:                      118.043
diff.sh:                       465.938  diff.sh:                       176.832  diff.sh:                       151.463
bi-grams.sh:                   752.060  bi-grams.sh:                   171.769  bi-grams.sh:                   156.463
set-diff.sh:                   1031.311 set-diff.sh:                   182.028  set-diff.sh:                   164.986
shortest-scripts.sh:           183.320  shortest-scripts.sh:           20.289   shortest-scripts.sh:           15.562
```

#### Unix50

#### Analytics MTS

```
1 336.147       1 84.240        1 49.946
2 336.735       2 85.818        2 50.433
3 408.517       3 94.363        3 60.014
4 107.901       4 68.830        4 35.527
```
