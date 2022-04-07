# OSDI'22 PaSh Artifact
>  This file is permanatly hosted as a [private gist (TODO)](XXXX). _Please do not share this file, as it includes account credentials to a server for running the experiments!_

## Questions/Notes/TODOs (for us -- delete before uploading):

- Do we want a video tutorial?
- Is the docker image pash:latest?
- Add graphviz in the quickcheck
- Shrink the current claim paragraphs and make them 1-2 sentences for each component and inline them in the list of components. KK: For now I have commented them out


This tutorial covers the `pash-jit`'s main functionality and key contributions as presented in the OSDI paper.
Here is the table of contents:

1. [Quick check](#quick-check)
2. [Components and Claims](#components-and-claims-artifact-functional)
3. TODO ...
<!-- 3. [Claim 1: Interfacing with the Shell](#claim-1-interfacing-with-the-shell)
4. [Claim 2: Stateful Parallelizing Compilation Server](#claim-2-a-stateful-parallelizing-compilation-server)
5. [Claim 3: Commutativity-aware Optimization](#claim-3-commutativity-aware-optimization)
6. [Experimental Evaluation](#experimental-evaluation) -->
7. [Support & Epilogue](#support--epilogue)

The recommended exploration path the following 
  (1) start with the [quick check](#quick-check), to confirm you can access everything (about 5 minutes)
  <!-- (2) watch the [video walkthrough](./artifact-video.mp4) of the entire artifact (about 20 minutes) -->
  (3) then cover details in this tutorial as necessary (the [full experimental run](#experimental-evaluation) takes several hours)

## Quick check

Here is how to check that everything works on your chosen evaluation environment:

#### Requirements

The artifact has the following requirements:
* Hardware: a modern multi-processor, to show performance results (the more cpus, the merrier)
* Software: automake bc curl gcc git libtool m4 python sudo wget bsdmainutils libffi-dev locales locales-all netcat-openbsd pkg-config python3 python3-pip python3-setuptools python3-testresources wamerican-insane

The `deathstar` server below meets all these requirements; 
  the `binpash/pash:latest` docker image is the next-best option, which meets all the software requirements.

#### Evaluation Server (deathstar)

We have created a reviewer account on `deathstar`, the machine used for all the performance results reported in the paper. Here _reviewers have to coordinate themselves other to not run checks/experiments at the same time_. To use `deathstar`, run a quick check with:

```sh
ssh osdi22@deathstar.ndr.md                         
./pash/scripts/quickcheck.sh                                           
```

#### Docker Container (local)

Another way to try PaSh is locally through a Docker container, running a pre-setup ubuntu Linux.
Information about docker installation may be found in [here](https://github.com/binpash/pash/tree/fixes/docs/install#docker-setup).
```sh
docker pull binpash/pash; docker run --name pash-playground -it binpash/pash
./opt/pash/scripts/quickcheck.sh                                            # this is typed _in_ the container
```


#### A Minimal Run: Demo Spell

All scripts in this guide assume that `$PASH_TOP` is set to the top directory of the PaSh codebase (i.e., `~/pash` on `deathstar` or `/opt/pash` in docker)

We will use `demo-spell.sh` --- a pipeline based [on the original Unix spell program](https://dl.acm.org/doi/10.1145/3532.315102) by Johnson --- to confirm that the infrastructure works as expected.

First, let's take a quick look at `spell`:

```sh
cd $PASH_TOP/evaluation/intro
cat demo-spell.sh
```

The script streams the input file into a pipeline that converts characters to lower case, removes punctuation, sorts in alphabetical order,  removes duplicate words, and filters out words from a dictionary file.

Next, let's run it on sequential inputs:

```sh
time ./demo-spell.sh > spell.out
```

We prefix the script with the `time` command, which should also output how long it took for the script to execute.
On `deathstar`, it takes about 41s.

To execute it using `pash` with 2x-parallelism:

```sh
time $PASH_TOP/pa.sh -w 2 -d 1 --log_file pash.log demo-spell.sh > pash-spell.out
``` 

On `deathstar`, the 2x-parallel script takes about 28s.

You can check that the results are correct by:

```sh
diff spell.out pash-spell.out
```

You could also execute it with 8x-parallelism using:
```sh
time $PASH_TOP/pa.sh -w 8 -d 1 --log_file pash.log demo-spell.sh > pash-spell.out
``` 

which takes about 14s.

To view the parallel code emitted by the compiler, you can inspect the log:

```sh
vim pash.log
```

The contents of the parallel script are shown after the line `(4) Will execute script in ...` and for 2x parallelism (`--width 2`) they should look like this:
```sh
rm_pash_fifos() {
...
mkfifo_pash_fifos() {
...
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ cat "input/100M.txt" >"/tmp/pash_khDRXr7/d8f7c89b966244cbb6b0d8dc0bfb9cb1/#fifo2" & }
pids_to_kill="${!} ${pids_to_kill}"
{ tr -cs A-Za-z "\\n" <"/tmp/pash_khDRXr7/d8f7c89b966244cbb6b0d8dc0bfb9cb1/#fifo4" >"/tmp/pash_khDRXr7/d8f7c89b966244cbb6b0d8dc0bfb9cb1/#fifo6" & }
pids_to_kill="${!} ${pids_to_kill}"
{ /home/osdi22/pash/runtime/auto-split.sh "/tmp/pash_khDRXr7/d8f7c89b966244cbb6b0d8dc0bfb9cb1/#fifo2" "/tmp/pash_khDRXr7/cee3744e05a64336a368295a133679e2/#fifo14" "/tmp/pash_khDRXr7/cee3744e05a64336a368295a133679e2/#fifo15" & }
pids_to_kill="${!} ${pids_to_kill}"
{ tr A-Z a-z <"/tmp/pash_khDRXr7/13d2518bd4ce45279a0128d64763b4f2/#fifo32" >"/tmp/pash_khDRXr7/cee3744e05a64336a368295a133679e2/#fifo17" & }
pids_to_kill="${!} ${pids_to_kill}"
{ tr A-Z a-z <"/tmp/pash_khDRXr7/cee3744e05a64336a368295a133679e2/#fifo15" >"/tmp/pash_khDRXr7/cee3744e05a64336a368295a133679e2/#fifo18" & }
pids_to_kill="${!} ${pids_to_kill}"
{ cat "/tmp/pash_khDRXr7/13d2518bd4ce45279a0128d64763b4f2/#fifo33" "/tmp/pash_khDRXr7/13d2518bd4ce45279a0128d64763b4f2/#fifo34" >"/tmp/pash_khDRXr7/d8f7c89b966244cbb6b0d8dc0bfb9cb1/#fifo4" & }
pids_to_kill="${!} ${pids_to_kill}"
{ /home/osdi22/pash/runtime/auto-split.sh "/tmp/pash_khDRXr7/d8f7c89b966244cbb6b0d8dc0bfb9cb1/#fifo6" "/tmp/pash_khDRXr7/cee3744e05a64336a368295a133679e2/#fifo19" "/tmp/pash_khDRXr7/cee3744e05a64336a368295a133679e2/#fifo20" & }
pids_to_kill="${!} ${pids_to_kill}"
{ sort <"/tmp/pash_khDRXr7/13d2518bd4ce45279a0128d64763b4f2/#fifo35" >"/tmp/pash_khDRXr7/cee3744e05a64336a368295a133679e2/#fifo22" & }
pids_to_kill="${!} ${pids_to_kill}"
{ sort <"/tmp/pash_khDRXr7/cee3744e05a64336a368295a133679e2/#fifo20" >"/tmp/pash_khDRXr7/cee3744e05a64336a368295a133679e2/#fifo23" & }
pids_to_kill="${!} ${pids_to_kill}"
{ sort -m "/tmp/pash_khDRXr7/13d2518bd4ce45279a0128d64763b4f2/#fifo36" "/tmp/pash_khDRXr7/13d2518bd4ce45279a0128d64763b4f2/#fifo37" >"/tmp/pash_khDRXr7/d8f7c89b966244cbb6b0d8dc0bfb9cb1/#fifo8" & }
pids_to_kill="${!} ${pids_to_kill}"
{ /home/osdi22/pash/runtime/auto-split.sh "/tmp/pash_khDRXr7/d8f7c89b966244cbb6b0d8dc0bfb9cb1/#fifo8" "/tmp/pash_khDRXr7/cee3744e05a64336a368295a133679e2/#fifo25" "/tmp/pash_khDRXr7/cee3744e05a64336a368295a133679e2/#fifo26" & }
pids_to_kill="${!} ${pids_to_kill}"
{ uniq <"/tmp/pash_khDRXr7/13d2518bd4ce45279a0128d64763b4f2/#fifo38" >"/tmp/pash_khDRXr7/cee3744e05a64336a368295a133679e2/#fifo28" & }
pids_to_kill="${!} ${pids_to_kill}"
{ uniq <"/tmp/pash_khDRXr7/cee3744e05a64336a368295a133679e2/#fifo26" >"/tmp/pash_khDRXr7/cee3744e05a64336a368295a133679e2/#fifo29" & }
pids_to_kill="${!} ${pids_to_kill}"
{ cat "/tmp/pash_khDRXr7/13d2518bd4ce45279a0128d64763b4f2/#fifo39" "/tmp/pash_khDRXr7/13d2518bd4ce45279a0128d64763b4f2/#fifo40" >"/tmp/pash_khDRXr7/cee3744e05a64336a368295a133679e2/#fifo30" & }
pids_to_kill="${!} ${pids_to_kill}"
{ uniq <"/tmp/pash_khDRXr7/cee3744e05a64336a368295a133679e2/#fifo30" >"/tmp/pash_khDRXr7/d8f7c89b966244cbb6b0d8dc0bfb9cb1/#fifo10" & }
pids_to_kill="${!} ${pids_to_kill}"
{ /home/osdi22/pash/runtime/eager.sh "/tmp/pash_khDRXr7/cee3744e05a64336a368295a133679e2/#fifo14" "/tmp/pash_khDRXr7/13d2518bd4ce45279a0128d64763b4f2/#fifo32" "/tmp/pash_khDRXr7/a0dfe82ede444443b9b3ec7de6dace82/pash_eager_intermediate_#fifo1" & }
pids_to_kill="${!} ${pids_to_kill}"
{ /home/osdi22/pash/runtime/eager.sh "/tmp/pash_khDRXr7/cee3744e05a64336a368295a133679e2/#fifo17" "/tmp/pash_khDRXr7/13d2518bd4ce45279a0128d64763b4f2/#fifo33" "/tmp/pash_khDRXr7/a0dfe82ede444443b9b3ec7de6dace82/pash_eager_intermediate_#fifo2" & }
pids_to_kill="${!} ${pids_to_kill}"
{ /home/osdi22/pash/runtime/eager.sh "/tmp/pash_khDRXr7/cee3744e05a64336a368295a133679e2/#fifo18" "/tmp/pash_khDRXr7/13d2518bd4ce45279a0128d64763b4f2/#fifo34" "/tmp/pash_khDRXr7/a0dfe82ede444443b9b3ec7de6dace82/pash_eager_intermediate_#fifo3" & }
pids_to_kill="${!} ${pids_to_kill}"
{ /home/osdi22/pash/runtime/eager.sh "/tmp/pash_khDRXr7/cee3744e05a64336a368295a133679e2/#fifo19" "/tmp/pash_khDRXr7/13d2518bd4ce45279a0128d64763b4f2/#fifo35" "/tmp/pash_khDRXr7/a0dfe82ede444443b9b3ec7de6dace82/pash_eager_intermediate_#fifo4" & }
pids_to_kill="${!} ${pids_to_kill}"
{ /home/osdi22/pash/runtime/eager.sh "/tmp/pash_khDRXr7/cee3744e05a64336a368295a133679e2/#fifo22" "/tmp/pash_khDRXr7/13d2518bd4ce45279a0128d64763b4f2/#fifo36" "/tmp/pash_khDRXr7/a0dfe82ede444443b9b3ec7de6dace82/pash_eager_intermediate_#fifo5" & }
pids_to_kill="${!} ${pids_to_kill}"
{ /home/osdi22/pash/runtime/eager.sh "/tmp/pash_khDRXr7/cee3744e05a64336a368295a133679e2/#fifo23" "/tmp/pash_khDRXr7/13d2518bd4ce45279a0128d64763b4f2/#fifo37" "/tmp/pash_khDRXr7/a0dfe82ede444443b9b3ec7de6dace82/pash_eager_intermediate_#fifo6" & }
pids_to_kill="${!} ${pids_to_kill}"
{ /home/osdi22/pash/runtime/eager.sh "/tmp/pash_khDRXr7/cee3744e05a64336a368295a133679e2/#fifo25" "/tmp/pash_khDRXr7/13d2518bd4ce45279a0128d64763b4f2/#fifo38" "/tmp/pash_khDRXr7/a0dfe82ede444443b9b3ec7de6dace82/pash_eager_intermediate_#fifo7" & }
pids_to_kill="${!} ${pids_to_kill}"
{ /home/osdi22/pash/runtime/eager.sh "/tmp/pash_khDRXr7/cee3744e05a64336a368295a133679e2/#fifo28" "/tmp/pash_khDRXr7/13d2518bd4ce45279a0128d64763b4f2/#fifo39" "/tmp/pash_khDRXr7/a0dfe82ede444443b9b3ec7de6dace82/pash_eager_intermediate_#fifo8" & }
pids_to_kill="${!} ${pids_to_kill}"
{ /home/osdi22/pash/runtime/eager.sh "/tmp/pash_khDRXr7/cee3744e05a64336a368295a133679e2/#fifo29" "/tmp/pash_khDRXr7/13d2518bd4ce45279a0128d64763b4f2/#fifo40" "/tmp/pash_khDRXr7/a0dfe82ede444443b9b3ec7de6dace82/pash_eager_intermediate_#fifo9" & }
pids_to_kill="${!} ${pids_to_kill}"
{ /home/osdi22/pash/runtime/eager.sh "/tmp/pash_khDRXr7/d8f7c89b966244cbb6b0d8dc0bfb9cb1/#fifo10" "/tmp/pash_khDRXr7/13d2518bd4ce45279a0128d64763b4f2/#fifo41" "/tmp/pash_khDRXr7/a0dfe82ede444443b9b3ec7de6dace82/pash_eager_intermediate_#fifo10" & }
pids_to_kill="${!} ${pids_to_kill}"
{ comm -13 input/sorted_words "/tmp/pash_khDRXr7/13d2518bd4ce45279a0128d64763b4f2/#fifo41" & }
pids_to_kill="${!} ${pids_to_kill}"
source /home/osdi22/pash/runtime/wait_for_output_and_sigpipe_rest.sh ${!} 2>>pash.log
rm_pash_fifos
...
```

Note that most stages in the pipeline are repeated twice and proceed in parallel (i.e., using `&`). This completes the "quick-check".

<!-- ## Video Tutorial -->


## Components and Claims (Artifact Functional)

Our paper describes the following system components:
<!-- kk: Here we should describe system components -->
* [Parsing library] ... TODO: Add link to parsing library
* [Preprocesor] -- The preprocessor [compiler/pash.py](https://github.com/binpash/pash/blob/main/compiler/pash.py) parses the shell script and instruments it with calls to [compiler/pash_runtime.sh](https://github.com/binpash/pash/blob/main/compiler/pash_runtime.sh).
* [JIT Engine] ... The JIT engine works through the interaction of the [compiler/pash_runtime.sh](https://github.com/binpash/pash/blob/main/compiler/pash_runtime.sh) instrumented calls and the stateful long lived compilation server. The runtime sends compilation requests to server and waits for a response:
  - If the server succeeds at compiling and parallizing the script, then the runtime runs the parallel shell script.
  - If the server fails, then it's not safe to parallize this region of the script and the runtime runs the original code.
* [Parallelizing Compilation Server] --  [compiler/pash_daemon.sh](https://github.com/binpash/pash/blob/main/compiler/pash_daemon.sh) handles compilation requests for parallelizing regions of the script. The compiler does hte following:
  - [Expansion] -- [compiler/expand.sh](https://github.com/binpash/pash/blob/main/compiler/expand.sh)
  - [Dependency untangling] -- If `--parallel_pipelines` flag is used, then the JIT Engine allows for running multiple regions at the same time which increases the amount of parallelization that can be extracted from the script. The safety check for running a new region is implemented [here](https://github.com/binpash/pash/blob/main/compiler/pash_runtime_daemon.py#:~:text=def-,check_resources_safety,-(self%2C%20process_id)%3A)
  - [Profile Driven] -- If `--profile-driven` flag is used, then the compiler adjusts the parallelization factor based on previous execution times. The responsible code for determining the configuration to use can be found [here](https://github.com/binpash/pash/blob/main/compiler/pash_runtime_daemon.py#:~:text=def-,determine_compiler_config,-(self%2C%20input_ir_file)%3A)
* [Commutativity Awareness] -- Consists of the following main components:
  - The annotations allow for adding the commutative property (e.g. [sort](https://github.com/binpash/pash/blob/main/annotations/sort.json))
  - The runtime nodes which are implemented in C for performance [c-split](https://github.com/binpash/pash/blob/main/runtime/r_split.c), [c-wrap](https://github.com/binpash/pash/blob/main/runtime/r_wrap.c), [c-strip](https://github.com/binpash/pash/blob/main/runtime/r_unwrap.c), and [c-merge](https://github.com/binpash/pash/blob/main/runtime/r_merge.c).
  - TODO: Should we have also the runtime.py commutativty changes to the transformation here?
<!-- * [Annotations](#claim-1-parallelizability-study--annotation-language): A study of the parallelizability of shell commands, and a lightweight annotation language for commands that are executable in a data-parallel manner.
* [Compiler](#claim-2-a-dataflow-based-parallelizing-compiler): A dataflow model and associated transformations that expose data parallelism while preserving the semantics of the sequential program.
* [Runtime](#claim-3-runtime-primitives): A set of runtime components, addressing the correctness and performance challenges arising during the execution of parallel shell scripts.

In the text below (but also more generally, in our discussions) we call first contribution "Annotations", the second contribution "Compiler", and the third contribution "Runtime". -->

And the following claims:
* [Verification of Dep Untangling] Add a link to the spin file and how to run it (TODO: KK)
* [Performance] ...
* [Correctness] ...

<!-- ## Claim 1: Interfacing with the Shell

The first claim is related to interface with the shell and is broken down into three parts. The first part is related to dynamic interposition and 
understanding the structure of ahead-of-time (AOT) parallization. The second part has to do with preprocessor and marking possible parallelizable regions 
with code that dynammicaly determine whether or not to invoke the compiler. The final part refers to PaSh-JIT parsing library during execution.

#### Dynamic Interposition

To understand the interposition techniques introduced by PaSh-JIT, we must first understand the simpler structure of ahead-of-time (AOT) parallelization.
The center of AOT parallelization is to indenty script fragments that may be safely parallelized to yield performance gains. The candidate fragments are called `pure` and should fullfill the following conditions:
1. Execute in parallel
2. Communicate with each other using explicit UNIX channels
3. Do not modify enviromental variables
4. They are fully expanded

PaSh-JIT works similar to AOT compiler and does the following steps:
1. Pause execution right before each parallelizable region
2. Compiling it to an efficient and queivalent parallel computation
3. Executing that instead

#### Preprocessor

Dynamic shell-script interposition without any shell-interpreter modifications is hard. To handle this, PaSh-JIT uses a light-weight script
instrumentation pre-processing step by marking candidate parallelizable regions and evaluation whether the invocation of the compiler is required.
However, PaSh-JIT analysis is impericese as there is no way to determine the pureness of an invoked command ahead of time. To address this,
PaSh-JIT does not explicitly search for parallelizable regions, but rather compilation sites. The preprocessor steps works by parsing lines of Bash
as they are consumed. The end goal of the preprocessor is to detect the maximal candidate parallelizable regions to:
1. minimize the number of compiler invocations
2. maximize the effect of each of them

#### Parsing Library

The final component related with shell interface, is the parsing/unparsing library. PaSh-JIT parses lines of shell-script as they are read and decomposes
the lines in order to execute them within the user's shell. Additinal calls to parsing/unparsing during compilation phase, when the server emits an 
optmized
Our parser is a reimplemented version of libdash in Python called pylibdash. This implementation uses Python bindings for the vanilla dash parser 
and reimplements the unparsing routines—adding a total of 0.9K LOC of Python over libdash.

## Claim 2: Stateful, Parallelizing Compilation Server


## Claim 3: Commutativity-aware Optimization -->

## Claim 1: Verification of Dependency Untangling

TODO: KK

## Claim 2: Correctness evaluation

TODO


## Claim 3: Performance evaluation



### Legend of terminology correspondence

Legend name correspondence between the paper and the artifact are seen below:
 - PaSh JIT               # OSDI: Blish
 - PaSh AOT               # OSDI: PaSh
 - PaSh JIT no_prof       # OSDI: Blish no_prof
 - PaSh JIT no_prof no_du # OSDI: Blish no_prof no_du
 - PaSh JIT no_comm       # OSDI: Blish no_comm

Flag name correspondence between the paper and the artifact are seen below:
 - PaSh JIT               # -w 16 --r_split --dgsh_tee --r_split_batch_size 1000000 --parallel_pipelines --profile_driven
 - PaSh AOT               # -w 16
 - PaSh JIT no_prof       # -w 16 --r_split --dgsh_tee --r_split_batch_size 1000000 --parallel_pipelines
 - PaSh JIT no_prof no_du # -w 16 --r_split --dgsh_tee --r_split_batch_size 1000000
 - PaSh JIT no_comm       # --parallel_pipelines --profile_driven



TODO: Maybe we could put that in a table?

- Flag names, legend names, (paper and code)

For example:
- dependency untangling -- parallel_pipelines

### Explain scripts, how much they are expected to take, and what do the figure that they produce look like (and why they are slightly different)

We offer two different input loads to test and evaluate the system; `--small` which corresponds to smaller inputs thus smaller execution time (4-5 hours) and `--full` input sizes that are used
in the paper (~20 hours; TODO we should verify that)

This section provides detailed instructions on how to replicate the figures of the experimental evaluation of the system as described in Table 1 from the paper (2-13).

![img](./imgs/table.png)


```sh
cd $PASH_TOP/evaluation/eval_script/
# There are two options here, either use --full as argument (default option) or --small (default option).
# This is the script that runs bash and PaSh JIT/AOT with several configuration. It runs all the benchmark described
# in Table 1 of the paper and gathers the execution times
bash run_all.sh --small
# after the execution is finished, you may view the raw data timers.
cat eval_results/run.tmp
# This will calculate the ratios of all the configurations in comparison with the bash execution timers (bash being the 
# baseline). It will produce the data_final.csv that will be used for generating the figures
bash gen_data.sh 
# this will produce figure{5,6,7}.pdf similar to the paper
./generate_charts.R data_final.csv
```

We have included in this repo sample data of the raw data timers (run.tmp), the final source data (data_final.csv) 
and the three output figures.

#### Pdfs
- [figure5](./pdfs/figure5.pdf)
- [figure6](./pdfs/figure6.pdf)
- [figure7](./pdfs/figure7.pdf)
 
#### Logs
 - [run.tmp](./logs/run.tmp)
 - [data_final.csv](./logs/data_final.csv)

#### Plot Generation

For the figures, additional packages are required (they are already available on `deathstar`):
```sh
# install the R environment and libraries
# follow the guide here https://cran.r-project.org/bin/linux/ubuntu/

# install ggplot2 for R
# needed for generating the pdf
sudo R
> install.packages('ggplot2', dep = TRUE)
> q()
> n
```

#### How to read Figures

- Explain legend, lines

Some results differences (with the paper):
- Some results are slightly worse w.r.t. absolute speedup (relative benefits are still the same) because the inputs are smaller and there is less potential for parallelizability
- In Figure 7, AvgTemp is faster in the artifact figure because we forgot to add parallel_pipelines (dependency untangling) and profile driven in the paper. Now we added them and therefore results are slightly better (relative speedups are still similar)

<!-- ## Experimental Evaluation

This section provides detailed instructions on how to replicate parts of the experimental evaluation of the system.

Note that input files that are used as inputs for this script are generated using the `gen*` scripts in [evaluation/scripts/input/](evaluation/scripts/input/).

```sh
# No need to run this on deathstar, inputs are already there
cd $PASH_TOP/evaluation/scripts/input/
./gen.sh
./gen_big_files.sh # Warning: This requires more than 100GB of space.
```

If you just want to run the scripts with small inputs (the main conclusions still hold)
you only need to run `./gen.sh`.

#### Section 6.1: Common Unix one-liners

The one-liner scripts are included in [evaluation/microbenchmarks](evaluation/microbenchmarks).
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

The script that runs PaSh on these programs is: [evaluation/eurosys/execute_eurosys_one_liners.sh](evaluation/eurosys/execute_eurosys_one_liners.sh) 
There are three modes of execution (can be seen by calling the script with the -h flag):

  1. Small inputs | --width 2, 16 | Only full PaSh config
  2. Small inputs | --width 2, 16 | All PaSh configs
  3. Big inputs | -- width 2, 4, 8, 16, 32, 64 | All PaSh configs

The script [evaluation/eurosys/execute_eurosys_one_liners.sh](evaluation/eurosys/execute_eurosys_one_liners.sh) is based on [evaluation/execute_compile_evaluation_script.sh](evaluation/execute_compile_evaluation_script.sh) that correctly sets up PaSh for the different configurations.

If you just want to check that PaSh achieves speedups as presented in the paper you can just run 1 with option `-s`.

If you are interested in seeing the improvements by PaSh's runtime primitives (all lines in Figure 9), you can run 2 with option `-m`. 
This should take a couple hours and should validate the trends between different PaSh configurations as shown in Figure 9.

If you want to reproduce the complete results from Figure 9, you need to run 3 with option `-l`.
Note that this should take more than a day to execute.
Also this requires several hundred GBs of free space (due to intermediate inputs, outputs, and buffering).

To plot the results from any of the above experiments, do the following:

```sh
cd $PASH_TOP/compiler
python3 gather_results.py
```

This will create plots for all invocations of [evaluation/eurosys/execute_eurosys_one_liners.sh`, one for each flag.
The plots are:
* for `-s`: [evaluation/plots/small_tiling_throughput_scaleup.pdf](evaluation/plots/small_tiling_throughput_scaleup.pdf)
* for `-m`: [evaluation/plots/medium_tiling_throughput_scaleup.pdf](evaluation/plots/medium_tiling_throughput_scaleup.pdf)
* for `-l`: [evaluation/plots/tiling_throughput_scaleup.pdf](evaluation/plots/tiling_throughput_scaleup.pdf)

Note that `-m` supersedes `-s` but `-l` does not supersede any of the two.

Also note that if you run a script partially, it might end up saving partial results,
therefore having 0 speedups in some points of the plots.

#### Section 6.2: Unix50 from Bell Labs

All of the Unix50 pipelines are in [evaluation/unix50/unix50.sh](evaluation/unix50/unix50.sh).
The inputs of the pipelines are in [evaluation/unix50/](evaluation/unix50/).

Before running the script we first need to move to the correct directory
  `cd $PASH_TOP/evaluation/eurosys`

The script that runs PaSh on these programs is: [evaluation/eurosys/execute_unix_benchmarks.sh](evaluation/eurosys/execute_unix_benchmarks.sh) 
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
python3 gather_results.py
```

This will create plots for both "1GB --width 4" and for "10GB --width 16".

The plots are in:
- for `-s`: [evaluation/plots/unix50_1GB_individual_speedups_4.pdf](evaluation/plots/unix50_1GB_individual_speedups_4.pdf)
- for `-l`: [evaluation/plots/unix50_10GB_individual_speedups_16.pdf](evaluation/plots/unix50_10GB_individual_speedups_16.pdf)

Note that the pipelines in the plot are sorted with respect to speedup, and not by their ID.
So the first pipeline does not necessarily correspond to the first pipeline in [evaluation/unix50](evaluation/unix50).

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

The program that we run, described in Section 6.3, can be seen in [evaluation/scripts/max-temp-complete.sh](evaluation/scripts/max-temp-complete.sh).
It takes as input a sequence of lines each containing a year (e.g. using `seq 2000 2004`).

To run the script with a single year of input use:
  `./execute_max_temp_dish_evaluation.sh -s`

These should take less than 10 minutes.

It runs the script on:
- bash
- pa.sh --width 16

The results are saved in:
- [evaluation/results/max-temp-complete-2000-2000-seq.time](evaluation/results/max-temp-complete-2000-2000-seq.time)
- [evaluation/results/max-temp-complete-2000-2000-16-pash.time](evaluation/results/max-temp-complete-2000-2000-16-pash.time)

If you want to run the program with 5 years of input (as is done in Section 6.3)
you need to use the following:
  `./execute_max_temp_dish_evaluation.sh -l`

It should take less than an hour. 
It also runs the script with bash and pash --width 16.

The results are saved in:
- [evaluation/results/max-temp-complete-2000-2004-seq.time](evaluation/results/max-temp-complete-2000-2004-seq.time)
- [evaluation/results/max-temp-complete-2000-2004-16-pash.time](evaluation/results/max-temp-complete-2000-2004-16-pash.time)

If you want to separate the preprocessing and processing (as done in Section 6.3)
you need to add the `-e` flag to either 1 or 5 year execution, e.g.:
  `./execute_max_temp_dish_evaluation.sh -l -e`

This runs:
- `evaluation/scripts/max-temp-preprocess.sh`
- `evaluation/scripts/max-temp-process.sh`

with bash, and pash --width 16. It saves results in:
- [evaluation/results/max-temp-preprocess-2000-2000-seq.time](evaluation/results/max-temp-preprocess-2000-2000-seq.time)
- [evaluation/results/max-temp-preprocess-2000-2000-16-pash.time](evaluation/results/max-temp-preprocess-2000-2000-16-pash.time)
- [evaluation/results/max-temp-process-2000-2000-seq.time](evaluation/results/max-temp-process-2000-2000-seq.time)
- [evaluation/results/max-temp-process-2000-2000-16-pash.time](evaluation/results/max-temp-process-2000-2000-16-pash.time)

and similarly for the large inputs (2000-2004).

Note that PaSh's speedup for the complete script 2000-2004 with width 16
is actually higher than what is reported in the paper since it doesn't
have to write the intermediate files (between preprocessing and processing) to disk.


#### Section 6.4: Use Case: Wikipedia Web Indexing

Note that input files that are needed by this script (complete Wikipedia) 
are saved locally on the server and therefore this program cannot be run from elsewhere.

Before running the script we first need to move to the correct directory
  `cd $PASH_TOP/evaluation/eurosys`

The program that we run, described in Section 6.4, can be seen in [evaluation/scripts/web-index.sh](evaluation/scripts/web-index.sh).
It requires having set the `$IN`, `$WIKI`, and `$WEB_INDEX_DIR` variables.

To run the script for a 1000 wikipedia links use:
  `./execute_web_index_dish_evaluation.sh -s`

This sets up the required variables and should take less than 5 minutes.
It runs the script with bash, pash --width 2, pash --width 16.

The results are saved in:
- [evaluation/results/web-index-1000-seq.time](evaluation/results/web-index-1000-seq.time)
- [evaluation/results/web-index-1000-2-pash.time](evaluation/results/web-index-1000-2-pash.time)
- [evaluation/results/web-index-1000-16-pash.time](evaluation/results/web-index-1000-16-pash.time)

If you want to run with the EuroSys evaluation inputs (100k links), use:
  `./execute_web_index_dish_evaluation.sh -l`

This should take a couple hours and the results are saved in:
- [evaluation/results/web-index-100000-seq.time](evaluation/results/web-index-100000-seq.time)
- [evaluation/results/web-index-100000-2-pash.time](evaluation/results/web-index-100000-2-pash.time)
- [evaluation/results/web-index-100000-16-pash.time](evaluation/results/web-index-100000-16-pash.time)

#### Section 6.5: Further Micro-benchmarks

To run the comparison with sort --parallel, just use [evaluation/eurosys/execute_baseline_sort.sh](evaluation/eurosys/execute_baseline_sort.sh)

Before running the script we first need to move to the correct directory
  `cd $PASH_TOP/evaluation/eurosys`

There are two modes of execution:
1. option: -s Small input | --width 2, 16
2. option: -l Big input | -- width 2, 4, 8, 16, 32, 64

Note that this script executes sort --parallel with double the value of --width
since we noticed that it grows slightly slower (as shown in the Figure in Section 6.5).

_This script throws a warning that is expected: `Env file: .../evaluation/microbenchmarks/sort_env_small.sh could not be read.` The warning is expected and can be safely ignored._ -->

## Support & Epilogue

PaSh is open-source and we are continuously working on it. Feel free to submit issues on Github, join the mailing lists, and use PaSh to parallelize your long-running shell scripts :)

### Useful Links:

**Website: TODO**

Mailing Lists: 
* [Discussion](https://groups.google.com/g/pash-discuss): Join this mailing list for discussing all things `pash`
* [Commits](https://groups.google.com/g/pash-commits): Join this mailing list for commit notifications

Development/contributions:
* Contribution guide: [docs/contrib](../../docs/contributing/contrib.md)
<!-- * Continuous Integration Server: [http://pash.ndr.md/](http://pash.ndr.md/) -->
