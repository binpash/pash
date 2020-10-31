# PaSh: Light-touch Data-Parallel Shell Processing

The repository contains:

* [compiler](./compiler): PaSh's data-flow model and associated optimization passes.
* [parallelizability](./parallelizability/): Parallelizability study used to inform the DSL characterizing commands.
* [evaluation](./evaluation): shell pipelines and example [scripts](./evaluation/scripts) used for the evaluation.
* [parser](./parser): The parser uses the `libdash` POSIX-compliant parser and `atdgen` to convert scripts to JSON-encoded ASTs.
* [paper](./paper): The paper submitted to EuroSys 2021, also available on [arxiv](https://arxiv.org/abs/2007.09436).
* [runtime](./runtime): A JVM-based runtime for PaSh (obsolete).
* [scripts](./scripts): Auxiliary scripts, for testing and running tasks.
* [tests](./tests): Contains a set of tests to check that everything works as expected (for CI).

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

# Tests
To execute the current tests (which are all the one-liners) simply run the following:

```sh
cd compiler
./test_evaluation_scripts.sh
```
