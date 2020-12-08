# PaSh: Light-touch Data-Parallel Shell Processing

The repository contains:

* [compiler](./compiler): PaSh's data-flow model and associated optimization passes.
* [parallelizability](./parallelizability/): Parallelizability study used to inform the DSL characterizing commands.
* [evaluation](./evaluation): shell pipelines and example [scripts](./evaluation/scripts) used for the evaluation.
* [parser](./parser): The parser uses the `libdash` POSIX-compliant parser and `atdgen` to convert scripts to JSON-encoded ASTs.
* [paper](./paper): The paper submitted to EuroSys 2021, also available on [arxiv](https://arxiv.org/abs/2007.09436).
* [runtime](./runtime): A JVM-based runtime for PaSh (obsolete).

## Installation

If you are on a fresh clone, first run the following to download the `libdash` submodule:
```sh
git submodule init
git submodule update
```

Then run:
```sh
./install.sh -p
```

Don't forget to export the library path in the end :)
```sh
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/"
```

## Tests
To execute the current tests (which are all the one-liners) simply run the following:

```sh
cd compiler
./test_evaluation_scripts.sh
```

## Running PaSh
To simply run PaSh on a script `script.sh` with parallelization width `2` make sure you are in the `compiler` directory and run:
```sh
python3.8 $PASH_TOP/compiler/pash.py --split_fan_out 2 script.sh
``` 

# Table of Contents

** Notes for artifact award: **
* create a video, basically outlining artifact evaluation and guiding EAC reviewers through all the results
* possibly create a minimal docker setup that people can test---with the right Python and Smoosh versions
* add README with full process and expected results
* optionally, split process into partial verification and full verification.

# High level Overview

PaSh is  a system for  parallelizing POSIX shell  scripts. Given a  script, PaSh
converts  it to  a dataflow  graph,  performs a  series of  semantics-preserving
program transformations that expose parallelism,  and then converts the dataflow
graph  back  into  a  script---one  that adds  POSIX  constructs  to  explicitly
guide parallelism  coupled with PaSh-provided Unix-aware  runtime primitives for
addressing performance- and correctness-related issues. A lightweight annotation
language allows  command developers to express  key parallelizability properties
about their commands.  An accompanying parallelizability study of  POSIX and GNU
commands—two large and  commonly used groups—guides the  annotation language and
optimized aggregator library that PaSh uses. PaSh’s extensive evaluation over 44
unmodified  Unix  scripts shows  significant  speedups  (0.89–61.1×, avg:  6.7×)
stemming  from  the  combination  of its  program  transformations  and  runtime
primitives.