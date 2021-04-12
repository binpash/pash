# A Short PaSh Tutorial

This tutorial covers the `pash`'s main functionality and key contributions as presented in the EuroSys paper.
Here is the table of contents:

0. [Introduction](#introduction)

1. [Installation](#installation)
2. [Running Introductory Scripts](#running-introductory-scripts)
3. [Experimental Evaluation](#experimental-evaluation)

4. [Parallelizability Study & Annotation Language](#parallelizability-study--annotation-language)
5. [A Dataflow-based Parallelizing Compiler](#a-dataflow-based-parallelizing-compiler)
6. [Runtime support](#runtime-support)

7. [Epilogue-development](#epilogue-development)

## Introduction

PaSh is a system for parallelizing POSIX shell scripts.
It has been shown to achieve order-of-magnitude performance improvements, while maintaining the correctness of the sequential scripts with respect to their sequential output.

> N.b.: PaSh is still under heavy development.

Consider the following log-analysis script, applied to many logs 

```sh
# log-analysis.sh
cat f1.md f2.md | 
  tr A-Z a-z |
  tr -cs A-Za-z '\n' |
  sort |
  uniq | 
  comm -13 dict.txt - > out
cat out | wc -l | sed 's/$/ mispelled words!/'
```
The first `cat` streams two markdown files into a pipeline that converts characters in the stream into lower case, removes punctuation, sorts the stream in alphabetical order, removes duplicate words, and filters out words from a dictionary file (line 1--7).
A second pipeline (line 7) counts the resulting lines to report the number of misspelled words to the user.

> If you're new to shell scripting, try to run each part of the pipeline separately and observe the output.
> For example, run `cat f1.md f2.md | tr A-Z a-z` on your terminal to witness the conversion to lower-case.

Visually, the script can be thought as executing as follows:

<img src="https://docs.google.com/drawings/d/e/2PACX-1vQv-Krzb9hxWCbbQC9Zg5knm6SySJrayh3mdZXG3Z4Y6hC4kgQj4PWqYmxNAR-LyKN5Fu9lWHJV0J0F/pub?w=517&amp;h=55">

The first pipeline (left; parts omitted) _sequentially_ processes `f1.md` and `f2.md` through all pipeline stages and writes to `out`.
After it, executes to completion the second pipeline starts its _sequential_ execution.

PaSh transforms and executes each pipeline in a data-parallel fashion.
Visually, the parallel script would look like this for 2x-parallelism (i.e., assuming that the computer on which we execute the script has at least two CPUs and that PaSh is invoked with `-w` value of `2`).

<img src="https://docs.google.com/drawings/d/e/2PACX-1vR8AL-gkL7CvqRJYiyX8z20_WcJ68l9JUEinJgI-_jKussb6q33Qlc1saaXx7Cf2pYp8-qjKhYXGu5e/pub?w=517&amp;h=55">

Given a script, PaSh converts it to a dataflow graph, performs a series of semantics-preserving program transformations that expose parallelism, and then converts the dataflow graph back into a POSIX script.
The new parallel script has POSIX constructs added to explicitly guide parallelism, coupled with PaSh-provided Unix-aware runtime primitives for addressing performance- and correctness-related issues.

#### PaSh Structure

PaSh consist of three main components and a few additional "auxiliary" files and directories. 
The three main components are:

* [annotations](../annotations/): DSL characterizing commands, parallelizability study, and associated annotations---more specifically, (i) a lightweight annotation language allows command developers to express key parallelizability properties about their commands; (ii) an accompanying parallelizability study of POSIX and GNU commands. guides the annotation language and optimized aggregator library 

* [compiler](../compiler): Shell-Dataflow translations and associated parallelization transformations---given a script, PaSh's compiler converts it to a dataflow graph, performs a series of semantics-preserving program transformations that expose parallelism, and then converts the dataflow graph back into a POSIX script. 

* [runtime](../runtime): Runtime components such as `eager`, `split`, and associated combiners. apart from POSIX constructs added to guide parallelism explicitly, PaSh provides Unix-aware runtime primitives for addressing performance- and correctness-related issues.

These three components implement the contributions presented [in the EuroSys paper](https://arxiv.org/pdf/2007.09436.pdf).
They are expected to be usable with minimal effort, through a few different installation means presented below.

The auxiliary directories are:
* [docs](../docs): Design documents, tutorials, installation instructions, etc.
* [evaluation](../evaluation): Shell pipelines and script used for the evaluation of `pash`.

## Installation

#### Ubuntu

If you're on ubuntu, setting up `pash` should be as easy as

```sh
curl up.pash.ndr.md | sh
```

#### Other Linux:

If you're on other linux distros, first install these packages (some of which might already be installed in your system):

```
libtool m4 automake pkg-config libffi-dev python3 python3-pip wamerican-insane bc bsdmainutils
```

Then clone the repo and run  `setup-pash.sh` as follows:
```sh
git clone git@github.com:andromeda/pash.git
./pash/script/setup-pash.sh
```

As noted at the end of `setup-pash.sh`, make sure you set `PASH_TOP` pointing to the absolute path of the directory `pa.sh` resides (you can optionally place that in your `PATH`).

#### Docker

Trying `pash` with Docker is easy.
Note that `pash` on Docker may or may not be able to exploit all available hardware resources, but is still useful for development on Windows or OS X.

**Installation:**
There are several options for installing `pash` via Docker.
The easiest is to `pull` the docker image [from GitHub](https://github.com/andromeda/pash/packages/715796):
```sh
docker pull docker.pkg.github.com/andromeda/pash/pash:latest
```
We refresh this image with every major release.

Alternatively, you can built the latest Docker container from scratch by running `docker build` in the repo:
```sh
git clone git@github.com:andromeda/pash.git
cd pash/scripts
docker build -t "pash/18.04" .
```
This will build a fresh Docker image using the latest commit---recommended for development.

It is also possible to fetch an image directly, but this is not longer recommended:
```sh
curl img.pash.ndr.md | docker load; docker run --name pash-play-$(whoami) -it pash/18.04
```

**Launch container**
After `pull`ing or `build`ing the image, you can launch the container as follows:
```sh
docker run --name pash-play -it pash/18.04
```

PaSh can be found in the container's `/pash` directory, so run `cd pash; git pull` to fetch the latest updates.
More information in the [pash-on-docker guide](./docs/contrib.md#pash-on-docker-a-pocket-guide).

## Running Scripts

All scripts in this guide assume that `$PASH_TOP` is set to the top directory of the PaSh codebase (e.g., `~/pash` on `deathstar` or `/pash` in docker)

> You can avoid including `$PASH_TOP` before `pa.sh` by adding `PASH_TOP` in your `PATH`, which amounts to adding an `export PATH=$PATH:$PASH_TOP` in your shell configuration file.

#### Hello World

The simplest script to try out `pash` is `hello-world`, which applies an expensive regular expression over the system dictionary file.

To run `hello-world.sh` normally, you would call `bash` on it:
```sh
time bash $PASH_TOP/evaluation/intro/hello-world.sh
```

To run it in parallel:
```sh
time $PASH_TOP/pa.sh $PASH_TOP/evaluation/intro/hello-world.sh
```

#### A More Interesting Script: Demo Spell

First, `cd` to `$PASH_TOP/evaluation/intro`:
```sh 
cd $PASH_TOP/evaluation/intro
```
We will use `demo-spell.sh` --- a pipeline based [on the original Unix spell program](https://dl.acm.org/doi/10.1145/3532.315102) by Johnson --- to confirm that the infrastructure works as expected. We need to setup the appropriate input files for this script to execute:
```sh
./input/setup.sh
```
After inputs are configured, let's take a quick look at `spell`:
```sh
cat demo-spell.sh
```
The script streams the input file into a pipeline that converts characters to lower case, removes punctuation, sorts in alphabetical order,  removes duplicate words, and filters out words from a dictionary file.

Next, let's run it on sequential inputs:
```sh
time ./demo-spell.sh > spell.out
```
We prefix the script with the `time` command, which should also output how long it took for the script to execute.
On our evaluation infrastructure, the script takes about 41s.

To execute it using `pash` with 2x-parallelism:
```sh
time $PASH_TOP/pa.sh -w 2 -d 1 --log_file pash.log demo-spell.sh > pash-spell.out
``` 
On `deathstar`, the 2x-parallel script takes about 28s.

You can check that the results are correct by:
```sh
diff spell.out pash-spell.out
```
Assuming you have more than 8 CPUs, you could also execute it with 8x-parallelism using:
```sh
time $PASH_TOP/pa.sh -w 8 -d 1 --log_file pash.log demo-spell.sh > pash-spell.out
``` 
which takes about 14s.

To view the parallel code emitted by the compiler, you can inspect the log:
```sh
cat pash.log
```

The contents of the parallel script are shown after the line `(4) Executing script in ...` and for 2x parallelism (`--width 2`) they should look like this:
```sh
rm -f "#file2"
...
mkfifo "#file2"
...
{ cat scripts/input/100M.txt >"#file2" & }
{ tr -cs A-Za-z "\\n" <"#file4" >"#file6" & }
{ /home/eurosys21/pash/runtime/auto-split.sh "#file2" "#file14" "#file15" & }
{ tr A-Z a-z <"#file32" >"#file17" & }
{ tr A-Z a-z <"#file15" >"#file18" & }
{ cat "#file33" "#file34" >"#file4" & }
{ /home/eurosys21/pash/runtime/auto-split.sh "#file6" "#file19" "#file20" & }
{ sort <"#file35" >"#file22" & }
{ sort <"#file20" >"#file23" & }
{ sort -m "#file36" "#file37" >"#file8" & }
{ /home/eurosys21/pash/runtime/auto-split.sh "#file8" "#file25" "#file26" & }
{ uniq <"#file38" >"#file28" & }
{ uniq <"#file26" >"#file29" & }
{ cat "#file39" "#file40" >"#file30" & }
{ uniq <"#file30" >"#file10" & }
{ /home/eurosys21/pash/runtime/eager.sh "#file14" "#file32" "/tmp/pash_eager_intermediate_#file1" & }
{ /home/eurosys21/pash/runtime/eager.sh "#file17" "#file33" "/tmp/pash_eager_intermediate_#file2" & }
{ /home/eurosys21/pash/runtime/eager.sh "#file18" "#file34" "/tmp/pash_eager_intermediate_#file3" & }
{ /home/eurosys21/pash/runtime/eager.sh "#file19" "#file35" "/tmp/pash_eager_intermediate_#file4" & }
{ /home/eurosys21/pash/runtime/eager.sh "#file22" "#file36" "/tmp/pash_eager_intermediate_#file5" & }
{ /home/eurosys21/pash/runtime/eager.sh "#file23" "#file37" "/tmp/pash_eager_intermediate_#file6" & }
{ /home/eurosys21/pash/runtime/eager.sh "#file25" "#file38" "/tmp/pash_eager_intermediate_#file7" & }
{ /home/eurosys21/pash/runtime/eager.sh "#file28" "#file39" "/tmp/pash_eager_intermediate_#file8" & }
{ /home/eurosys21/pash/runtime/eager.sh "#file29" "#file40" "/tmp/pash_eager_intermediate_#file9" & }
{ /home/eurosys21/pash/runtime/eager.sh "#file10" "#file41" "/tmp/pash_eager_intermediate_#file10" & }
{ comm -13 scripts/input/dict.txt "#file41" & }
source /home/eurosys21/pash/runtime/wait_for_output_and_sigpipe_rest.sh ${!}
rm -f "#file2"
...
```

Note that most stages in the pipeline are repeated twice and proceed in parallel (i.e., using `&`). This completes the "quick-check".

## Parallelizability Study & Annotation Language

PaSh includes (i) a parallelizability study of commands in POSIX and GNU Coreutils, and (ii) an annotation language for describing the parallelizability properties of individual commands.
The parallelizability study informed the design of the annotation language, which was in turn used to capture the key parallelizability characteristics in many of these commands.

#### Parallelizability Study of Commands in GNU & POSIX

The main results of the parallelizability study are summarized in the paper (Sec. 3.1 and Tab. 1).
In the artifact, the parallelizability study is summarized in the [./annotations/p_stats](../annotations/p_stats).

We note that, like the rest of the project, the results of the parallelizability study are being updated and might not exactly reflect percentages in the table;
  the camera-ready version of the paper will include the final numbers.

#### Annotations of Commands in GNU & POSIX

Annotations can be thought of as defining a bidirectional correspondence between a command and a node in the dataflow graph.
Since command behaviors (and correspondence) can change based on their arguments, annotations contain a sequence of predicates.
Each predicate is accompanied by information that instantiates the correspondence between a command and a dataflow node.

Annotations for about 60 popular commands are stored in [./annotations](../annotations) encoded as JSON.
These average about 14 lines per annotation, for a total of 846 lines of annotations.

Below we present two example annotations for `chmod` and `cut`.

```json
{
  "command": "chmod",
  "cases": [
    {
      "predicate": "default",
      "class": "side-effectful"
    }
  ]
}
```

The annotation for `chmod` is very simple, since it only needs to establish that `chmod` is side-effectful and therefore cannot be translated to a dataflow node.

```json
{
  "command": "cut",
  "cases": [
    {
      "predicate": {
        "operator": "or",
        "operands": [
          {
            "operator": "val_opt_eq",
            "operands": [
              "-d",
              "\n"
            ]
          },
          {
            "operator": "exists",
            "operands": [
              "-z"
            ]
          }
        ]
      },
      "class": "pure",
      "inputs": [
        "args[:]"
      ],
      "outputs": [
        "stdout"
      ]
    },
    {
      "predicate": "default",
      "class": "stateless",
      "inputs": [
        "args[:]"
      ],
      "outputs": [
        "stdout"
      ]
    }
  ],
  "options": [
    "stdin-hyphen",
    "empty-args-stdin"
  ],
  "short-long": [
    {
      "short": "-d",
      "long": "--delimiter"
    },
    {
      "short": "-z",
      "long": "--zero-terminated"
    }
  ]
}
```

The annotation for `cut` has two cases, each of which consists of a predicate on its arguments, and then an assignment of its parallelizability class, inputs, and outputs.
The first predicate indicates that cut is “pure”, i.e., not parallelizable but representable as a dataflow node, if the value accompanying the `-d` option is `\n` or if it was used with the `-z` flag.
In both of these cases, newlines do not represent data item boundaries, but are rather used internally by the command, making it unsafe to parallelize by splitting on line boundaries.
In all other cases (see the “default” case) the command is stateless.
Inputs are always assigned to the non-option arguments and the output is always stdout.
The option “stdin-hyphen” indicates that a non-option argument that is just a dash `-` represents the stdin, and the option “empty-args-stdin” indicates that if non-option arguments are empty, then the command reads from its stdin.
The list identified by “short-long” contains a correspondence of short and long argument names for this command.

Note that currently PaSh's `cut` annotation engine does not handle the `empty-args-stdin` option and the `short-long` key, so the annotation looks slightly different. We are continuously working on the annotations engine to improve it to handle a wider range of annotations. 

## A Dataflow-based Parallelizing Compiler

The bulk of PaSh is in the [compiler](../compiler/) directory.
A diagram of the compiler is shown in [docs/pash_architecture.png](../docs/pash_architecture.png).
A correspondence between blocks in the diagram and Python modules is shown below:

- PaSh Preprocessor -- [pash.py](../compiler/pash.py)
- Expand, Compile -- [ast_to_ir.py](../compiler/ast_to_ir.py)
- Annotations -- [annotations.py](../compiler/annotations.py), [command_categories.py](../compiler/command_categories.py)
- Optimize -- [pash_runtime.py](../compiler/pash_runtime.py)

**Note:** At the time of the paper submission, PaSh did not have a preprocessing component, and didn't handle variable expansion. These changes significantly improve the practical applicability of PaSh since it can be used on scripts where the environment variables are modified throughout the script.

First, there is the parser in [compiler/parser](../compiler/parser), which is a port of [libdash](https://github.com/mgree/), the dash parser extended with OCaml bindings, extended with ocaml2json and json2ocaml code to interface with PaSh.

Now let's get to the compiler. It's entry point is [compiler/pash.py](../compiler/pash.py) that parses a script and replaces potentially parallelizable regions with calls to [compiler/pash_runtime.sh](../compiler/pash_runtime.sh). It then executes the script.
This allows invoking the compiler during the runtime to have information about the values of environment variables.

The runtime script [compiler/pash_runtime.sh](../compiler/pash_runtime.sh) simply invokes the compiler [compiler/pash_runtime.py](../compiler/pash_runtime.py) and if it succeeds it executes the optimized script, otherwise it executes the original script.

Now the compiler has several stages:

1. It expands words in the AST and then it turns it into our dataflow model (guided by annotations)
   - The expansion and translation happens in [ast_to_ir.py](../compiler/ast_to_ir.py)
   - The dataflow model is mostly defined in [ir.py](../compiler/ir.py)
   - The annotations are processed in [annotations.py](../compiler/annotations.py) and [command_categories.py](../compiler/command_categories.py)
2. It performs transformations on the dataflow graph to expose parallelism (guided by annotations)
   - Translations happen in [pash_runtime.py](../compiler/pash_runtime.py)
3. It then translates the dataflow graph back to a shell script to execute it with bash
   - The `dfg2shell` translation happens in [ir_to_ast.py](../compiler/ir_to_ast.py)
   
 A few interesting fragments are shown below.
 
 The [ast_to_ir.py](https://github.com/andromeda/pash/blob/main/compiler/ast_to_ir.py) contains a case statement that essentially pattern-matches on constructs of the shells script AST and then compiles them accordingly.
```Python
 compile_cases = {
        "Pipe": (lambda fileIdGen, config:
                 lambda ast_node: compile_node_pipe(ast_node, fileIdGen, config)),
        "Command": (lambda fileIdGen, config:
                    lambda ast_node: compile_node_command(ast_node, fileIdGen, config)),
        "And": (lambda fileIdGen, config:
                lambda ast_node: compile_node_and_or_semi(ast_node, fileIdGen, config)),
        "Or": (lambda fileIdGen, config:
               lambda ast_node: compile_node_and_or_semi(ast_node, fileIdGen, config))
        # ... more code ...
    
```


The following function from [ir.py](https://github.com/andromeda/pash/blob/main/compiler/ir.py) is responsible for parallelizing a single node (i.e., command) in the dataflow graph. Look at the schematic in the comments starting [on line 637](https://github.com/andromeda/pash/blob/main/compiler/ir.py#L637) that gives the high-level overview of what this function does (not shown below).

```Python
    # See comment on line 637
    def parallelize_node(self, node_id, fileIdGen):
        node = self.get_node(node_id)
        assert(node.is_parallelizable())

        ## Initialize the new_node list
        new_nodes = []
        # ... more code ...
```

Another interesting fragment is in [ir_to_ast.py](https://github.com/andromeda/pash/blob/main/compiler/ir_to_ast.py), which translates the parallel dataflow graph back to an AST.

```Python
def ir2ast(ir, args):
    # ... more code ...
    
    ## Make the main body
    body = ir.to_ast(drain_streams)
    
    # ... more code ...
    
    ## Call the prologue that creates fifos for all ephemeral fids    
    prologue = make_ir_prologue(ephemeral_fids)
    
    ## Call the epilogue that removes all ephemeral fids
    epilogue = make_ir_epilogue(ephemeral_fids, clean_up_graph, args.log_file)

    final_asts = prologue + body + epilogue
```

This AST is then unparsed back into a (parallel) shell script.

## Runtime primitives

PaSh includes a small library of runtime primitives for supporting the runtime execution of parallel scripts emitted by the compiler.

### Stream Splitting

The PaSh compiler inserts `split` nodes to expose parallelism when parallelizable nodes only have one input.
The command `split` is in [./runtime](../runtime).

### Eager Stream Polling

To overcome the laziness challenges outlined in Sec. 5, PaSh inserts and instantiates `eager` nodes on streams.
The command `eager` is in [./runtime](../runtime).

### Cleanup Logic

PaSh contains cleanup logic for dealing with dangling FIFOs.
This is implemented in `wait_for_output_and_sigpipe_rest.sh`, contained in [./runtime](../runtime).

### Aggregators

There is a small custom aggregator library provided in [runtime/agg/py/](../runtime/agg/py/).
These aggregators are used to merge partial results from the parallel scripts.

For example, the aggregator `wc.py` can merge results from partial `wc` running in parallel.
To confirm what the aggregator does, call it as follows.

```shell 
$PASH_TOP/runtime/agg/py/wc.py <(echo -e '1\n2\n3' | wc -l) <(echo -e '1\n2\n3' | wc -l)
```

Internally, the aggregator looks like this:

```Python
#!/usr/bin/python
import sys, os, functools, utils

def parseLine(s):
  return map(int, s.split())

def emitLine(t):
  return [" ".join(map(lambda e: str(e).rjust(utils.PAD_LEN, ' '), t))]

def combiner(a, b):
  if not a:
    return b
  az = parseLine(a[0])
  bz = parseLine(b[0])
  return emitLine([ (i+j) for (i,j) in zip(az, bz) ])

utils.help()
res = functools.reduce(agg, utils.read_all(), [])
utils.out("".join(res))
```

The core of the aggregator, function `combiner`, is binary (i.e., takes two input streams).
The `reduce` function lifts the aggregator to arity _n_, i.e., the number of the incoming edges—each of which feeds the aggregator with the results of running a parallel `wc` instance.
This allows developers to think of aggregators in terms of two inputs, but generalize their aggregators to many inputs.
Utility functions such as `read` and `help`, common across our aggregator library, deal with error handling when reading multiple file descriptors, and offer a invocation flag `-h` that demonstrates the use of each aggregator.

PaSh’s library currently several aggregators, many of which are usable by more than one command or flag. For example, the aggregator shown above is shared among `wc`, `wc -lw`, `wc -lm` etc.

## Experimental Evaluation

_Most benchmark sets in the evaluation infrastructure include a `input/setup.sh` script for fetching inputs and setting up the experiment appropriately._
See [Running other script]() later.

This section provides detailed instructions on how to replicate parts of the experimental evaluation of the system.

Note that input files that are used as inputs for this script are generated using the `gen*` scripts in [evaluation/scripts/input/](../evaluation/scripts/input/).

```sh
# No need to run this on deathstar, inputs are already there
cd $PASH_TOP/evaluation/scripts/input/
./gen.sh
./gen_big_files.sh # Warning: This requires more than 100GB of space.
```

If you just want to run the scripts with small inputs (the main conclusions still hold)
you only need to run `./gen.sh`.

#### Section 6.1: Common Unix one-liners

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

#### Section 6.2: Unix50 from Bell Labs

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

#### Section 6.3: Use Case: NOAA Weather Analysis

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


#### Section 6.4: Use Case: Wikipedia Web Indexing

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


## Support & Epilogue

PaSh is open-source and we are continuously working on it, hoping that it will "escape" the research prototype orbit and be useful in practice. Feel free to submit issues on Github, join the mailing lists, and use PaSh to parallelize your long-running shell scripts :)

### Useful Links:

Mailing Lists: 
* [Discussion](https://groups.google.com/g/pash-discuss): Join this mailing list for discussing all things `pash`
* [Commits](https://groups.google.com/g/pash-commits): Join this mailing list for commit notifications

Development/contributions:
* Contribution guide: [docs/contrib](../docs/contrib.md)
* Continuous Integration Server: [http://pash.ndr.md/](http://pash.ndr.md/)
