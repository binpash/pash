# OSDI'22 PaSh Artifact
>  This file is permanently hosted as a [private gist (TODO)](XXXX). _Please do not share this file, as it includes account credentials to a server for running the experiments!_

## TODOs (for us -- delete before uploading):

- Move this to a private gist and add the keys.

- Is the docker image pash:latest? No, Docker images update automatically on each new release. I could try updating them manually

- Make a new branch with these fixes (called osdi22) and freeze it. We can then merge this branch to `future`. Once we do that, we also need to fix links to the `fixes` branch to point to the new branch.

- In general check all links once we are done. This will be a gist so relative links won't work.

- **Thurston:** Check that the link/description of the parsing library are OK.

This artifact covers the `PaSh-JIT`'s main functionality and key contributions as presented in the OSDI paper.
Here is the table of contents of this README file:

1. [Quick check](#quick-check)
2. [System Components](#components-artifact-functional)
3. [Claim 1 -- Correctness Evaluation](#claim-1-correctness-evaluation)
4. [Claim 2 -- Performance Evaluation](#claim-2-performance-evaluation)
5. [Claim 3 -- Dependency Untangling Verification](#claim-3-verification-of-dependency-untangling)
6. [Support & Epilogue](#support--epilogue)

The recommended artifact evaluation path is the following:

1. Start with the [quick check](#quick-check), to confirm you can access everything (less than 20 minutes).

2. Read a brief overview of the [system components](#components-artifact-functional) that includes pointers to the code that implements them.

3. Evaluate the three evaluation claims that are made in the paper, namely 
  (i) the [correctness evaluation](#claim-1-correctness-evaluation), 
    which corresponds to Paper Section 7.1,
  (ii) the [performance evaluation](#claim-2-performance-evaluation) (takes multiple hours), 
    which corresponds to Paper Section 7.2, and
  (iii) the [verificaton of the dependency untangling algorithm](#claim-3-verification-of-dependency-untangling), which is mentioned in the end of Section 5.2.

## Quick check

Here is how to check that everything works on your chosen evaluation environment:

### Installation/Infrastructure

#### Requirements

The artifact makes three claims, each of which has different hardware/software requirements:
- Claim 1 requires access to one of our machines in order to run the POSIX test suite because it cannot be shared publicly.
- Claim 2 has the following requirements, which are met by a server that we give you access to, or can be met using the `binpash/pash:latest` docker image (but you would have to run it on your own infrastructure):
  * Hardware: a modern multi-processor, to show performance results (the more cpus, the merrier)
  * Software: automake bc curl gcc git libtool m4 python sudo wget bsdmainutils libffi-dev locales locales-all netcat-openbsd pkg-config python3 python3-pip python3-setuptools python3-testresources wamerican-insane
- Claim 3 requires **TODO**


#### Correctness Evaluation Server (antikythera)

**TODO: Dimitris** Have a section that runs quickcheck (or something small like `./pa.sh -c 'echo hi'`) here (together with instructions for keys).

#### Performance Evaluation Server (deathstar)

We have created a reviewer account on `deathstar`, the machine used for all the performance results reported in the paper. Here _reviewers have to coordinate themselves other to not run checks/experiments at the same time_. To use `deathstar`, run a quick check with:

```sh
ssh osdi22@deathstar.ndr.md
# this will print information regarding pash installation
./pash/scripts/quickcheck.sh                                           
```

**TODO:** What information will this print? Can we paste some (all) of it here?

##### Private Key Installation

**TODO: Dimitris** add instructions on how to install the private ssh key



#### Docker Container (local)

Another way to try PaSh is locally through a Docker container, running a pre-setup ubuntu Linux.
Information about docker installation may be found in [here](https://github.com/binpash/pash/tree/fixes/docs/install#docker-setup).
```sh
docker pull binpash/pash; docker run --name pash-playground -it binpash/pash
./opt/pash/scripts/quickcheck.sh                                            # this is typed _in_ the container
```

##### Plot Generation

In order to generate the performance evaluation figures, additional packages are required (they are already available on `deathstar` but not on the docker image):
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


### A Minimal Run: Demo Spell

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
# download the input data
cd input && bash setup.sh && cd ../
time ./demo-spell.sh > spell.out
```

We prefix the script with the `time` command, which should also output how long it took for the script to execute.
On `deathstar`, it takes about 41s.

To execute it using `pash` with 2x-parallelism:

```sh
time $PASH_TOP/pa.sh -w 2 demo-spell.sh > pash-spell.out
``` 

On `deathstar`, the 2x-parallel script takes about 28s.


You can check that the results are correct by:

```sh
cmp spell.out pash-spell.out && echo "Files are identical"
```

You could also execute it with 8x-parallelism using:
```sh
time $PASH_TOP/pa.sh -w 8 --graphviz pdf --graphviz_dir . demo-spell.sh > pash-spell.out
``` 
which takes about 14s.

The `--graphviz` option generates a graph in a directory `./pash_graphviz_<TIMESTAMP>`.

A sample graph may be viewed [here](./pdfs/dfg.pdf).



To view the actual parallel shell code emitted by the compiler, you can inspect the log by running:

```sh
time $PASH_TOP/pa.sh -w 8 -d 1 --log_file pash.log demo-spell.sh > pash-spell.out
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

...

{ comm -13 input/sorted_words "/tmp/pash_khDRXr7/13d2518bd4ce45279a0128d64763b4f2/#fifo41" & }
pids_to_kill="${!} ${pids_to_kill}"
source /home/osdi22/pash/runtime/wait_for_output_and_sigpipe_rest.sh ${!} 2>>pash.log
rm_pash_fifos
...
```

Note that most stages in the pipeline are repeated twice and proceed in parallel (i.e., using `&`). This completes the "quick-check".


## Components (Artifact Functional)

Our paper describes the following system components. An overview of how they all interact can be found in the Paper Figure 1.
* **Parsing library:** The parsing library (Paper Section 3.3) contains Python bindings for the dash parser (e.g., [ast2a.py](https://github.com/binpash/pash/blob/main/compiler/parser/ceda/ast2a.py) translates dash's AST to a Python AST) and a complete [unparser implementation](https://github.com/binpash/pash/blob/main/compiler/parser/ceda/ast2shell.py).
* **Preprocesor:** The preprocessor [compiler/pash.py](https://github.com/binpash/pash/blob/main/compiler/pash.py) (Paper Section 3.2) uses the parser and then instruments the AST with calls to the [JIT Engine](https://github.com/binpash/pash/blob/main/compiler/pash_runtime.sh).
* **JIT Engine:** The JIT engine (Paper Section 4) [compiler/pash_runtime.sh](https://github.com/binpash/pash/blob/main/compiler/pash_runtime.sh) transitions between shell and PaSh mode and interacts with the the stateful compilation server. The engine sends compilation requests to server and waits for a response:
  - If the server succeeds at compiling and parallelizing the requested region, then the engine runs the parallel shell script.
  - If the server fails, then it's not safe to parallelize this region and the engine runs the original code.
* **Parallelizing Compilation Server:** The parallelizing compilation server (Paper Section 5) [compiler/pash_runtime_daemon.py](https://github.com/binpash/pash/blob/main/compiler/pash_runtime_daemon.py) handles compilation requests for parallelizing regions of the script. The server contains the following subcomponents:
  - **Expansion:** The expansion module (Paper Section 5.1) that expands script fragments [compiler/expand.py](https://github.com/binpash/pash/blob/main/compiler/expand.py).
  - **Dependency untangling:** The dependency untangling component (Paper Section 5.2) is invoked using the `--parallel_pipelines` flag. When enabled, the JIT Engine allows for running multiple regions at the same time; This increases the amount of parallelization that can be extracted from the script. The safety check for detecting dependencies between fragments [pash_runtime_daemon.py, L155](https://github.com/binpash/pash/blob/ebe379427a42040c1c00b01bcdcadb1fdbfb0281/compiler/pash_runtime_daemon.py#L155) is called after compilation succeeds [pash_runtime_daemon.py, L317](https://github.com/binpash/pash/blob/d836bfd32a58c848dccc2b157e72237f412c6ff5/compiler/pash_runtime_daemon.py#L317). Finally, the code responsible for executing the optimized script in parallel is in the JIT Engine [pash_runtime.sh, L257](https://github.com/binpash/pash/blob/d836bfd32a58c848dccc2b157e72237f412c6ff5/compiler/pash_runtime.sh#L257)
  - **Profile Driven Configuration:** The profile driven component (Paper Section 5.3) is invoked using the `--profile-driven`. When enabled, the compiler adjusts the parallelization factor based on previous execution times. The code that determines the parallelization based on profiles is in [pash_runtime_daemon.py, L172](https://github.com/binpash/pash/blob/ebe379427a42040c1c00b01bcdcadb1fdbfb0281/compiler/pash_runtime_daemon.py#L172).
* **Commutativity Awareness:** The commutativity-aware components of PaSh-JIT (Paper Section 6) consist of the following main components:
  - The annotations allow for indicating that a command is commutative (e.g., [sort](https://github.com/binpash/pash/blob/ebe379427a42040c1c00b01bcdcadb1fdbfb0281/annotations/sort.json#L18))
  - The dataflow nodes for orchestrating commutativity-aware concurrency are [c-split](https://github.com/binpash/pash/blob/main/runtime/r_split.c), [c-wrap](https://github.com/binpash/pash/blob/main/runtime/r_wrap.c), [c-strip](https://github.com/binpash/pash/blob/main/runtime/r_unwrap.c), and [c-merge](https://github.com/binpash/pash/blob/main/runtime/r_merge.c).


## Claim 1: Correctness evaluation

This corresponds to Section 7.1 of the paper.

### Requirements

In order to run the correctness evaluation, you need to login the `antikythera` machine since the POSIX test suite cannot be shared publicly.

We have already installed PaSh and the POSIX suite in a docker image there that you can use to run the tests and validate the correctness results.

**TODO: Dimitris** Add command/instructions to ssh into antikythera with the key `-i ...`.

```sh
# ssh to antikythera.csail.mit.edu using the provided private key
# inside the machine, run
docker run -it --sig-proxy=false posix /bin/bash
# this will spawn a clean docker instance to run the posix tests using bash and pash
# you should then be greeted with this message:
# osdi22@antikythera:~$ docker run -it --sig-proxy=false posix /bin/bash
# using TET_ROOT = /home/runner/tet3.8
# VSC environment setup is successful

# First, setup the posix tests to pash configuration
bash ~/prepare_pash.sh
cd ~/tet3.8/vsc
# Prepare the tests
tcc -bp vsc posix_shell
# Run the tests (~10 minutes)
tcc -ep vsc posix_shell

# You can use the following script to view the detailed results
bash ~/results/summarize_journal.sh ~/tet3.8/vsc/results/0002e/journal

## Which should finish like this:

#  ...
#  ...
# ========================================================================
#  494 tests:
#  375         passed
#   41         failed
#   31       untested
#    6     unresolved
#   40    unsupported
#    1     not in use
#    0   other status

# Then, setup the posix tests to use bash
bash ~/prepare_bash.sh
cd ~/tet3.8/vsc
# Prepare the tests
tcc -bp vsc posix_shell
# Run the tests (~10 minutes)
tcc -ep vsc posix_shell

# You can use the following script to view the generated bash results
bash ~/results/summarize_journal.sh ~/tet3.8/vsc/results/0004e/journal

## Which should finish like this:

#  ...
#  ...
# ========================================================================
#  494 tests:
#  377         passed
#   39         failed
#   31       untested
#    6     unresolved
#   40    unsupported
#    1     not in use
#    0   other status

## You can then run the following command to see the difference of two tests:
comm -13 <(../../results/summarize_journal.sh results/0002e/journal | grep "Assertion #" | sort) <(../../results/summarize_journal.sh results/0004e/journal | grep "Assertion
#" | sort)
## This should return the following
# Assertion #430 (A): When a command fails during word expansion or redirection, then
# Assertion #691 (A): When the shell is not executing interactively, then the 'set -u'

## which are exactly the tests mentioned in 
```

## Claim 2: Performance evaluation

This corresponds to the section 7.2 of the paper.

### Legend of terminology correspondence

Flag name correspondence between the paper and the artifact are seen below:
 - PaSh JIT               # -w 16 --r_split --dgsh_tee --r_split_batch_size 1000000 --parallel_pipelines --profile_driven
 - PaSh AOT               # -w 16
 - PaSh JIT no_prof       # -w 16 --r_split --dgsh_tee --r_split_batch_size 1000000 --parallel_pipelines
 - PaSh JIT no_prof no_du # -w 16 --r_split --dgsh_tee --r_split_batch_size 1000000
 - PaSh JIT no_comm       # --parallel_pipelines --profile_driven


### Executing the evaluation and plotting the results

We offer two different input loads to test and evaluate the system; `--small` which corresponds to smaller inputs thus smaller execution time (4-5 hours) and `--full` input sizes that are used
in the paper (>20 hours). Running the small results returns results that are very close to the ones shown in the paper, and all differences between configurations are evident. The only difference is that a few speedups are slightly smaller (for scripts that are too small to get meaningful benefits from parallelization).

This section provides detailed instructions on how to replicate the figures of the experimental evaluation of the system as described in Table 1 from the paper (2-13).

![img](./imgs/table.png)

**TODO: Dimitris** Point to each script from the table to a source file.

```sh
cd $PASH_TOP/evaluation/eval_script/
# There are two options here, either use --small or --full as an argument to determine the input size.
# This is the script that runs bash and PaSh JIT/AOT with several configuration. It runs all the benchmark described
# in Table 1 of the paper and gathers the execution times

# We suggest that you use `tmux/screen` to run this script so that you can leave it running in the background and see the results after it is done.
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
- [Paper, Figure 5](./pdfs/figure5.pdf)
- [Paper, Figure 6](./pdfs/figure6.pdf)
- [Paper, Figure 7](./pdfs/figure7.pdf)
 
#### Logs
 - [run.tmp](./logs/run.tmp)
 - [data_final.csv](./logs/data_final.csv)


#### Interpreting the figures

The figures are slightly different from the ones shown in the paper for the following two reasons:
- For some benchmarks, the smaller input leads to a more pronounced cost of one-time overheads and therefore smaller absolute speedup for all PaSh configurations (JIT, AOT, etc). However, the relative benefits between configurations are still the same.
- In Figure 7, AvgTemp is faster in the artifact figure because in the paper we forgot to use `--parallel_pipelines` (dependency untangling) and `--profile_driven` when we run the experiment for the paper. Now we added them and therefore results are slightly better (relative speedups are still similar).

## Claim 3: Verification of Dependency Untangling

TODO: KK


## Support & Epilogue

PaSh is open-source and we are continuously working on it. Feel free to submit issues on Github, join the mailing lists, and use PaSh to parallelize your long-running shell scripts :)

### Useful Links:

Website: [binpa.sh](https://binpa.sh)

Mailing Lists: 
* [Discussion](https://groups.google.com/g/pash-discuss): Join this mailing list for discussing all things `pash`
* [Commits](https://groups.google.com/g/pash-commits): Join this mailing list for commit notifications

Development/contributions:
* Contribution guide: [docs/contrib](../../docs/contributing/contrib.md)
<!-- * Continuous Integration Server: [http://pash.ndr.md/](http://pash.ndr.md/) -->
