# A Short PaSh Tutorial
Quick jump: [Introduction](#introduction) | [Running Scripts](#running-scripts) | [What Next?](#what-next)

This short tutorial covers the `pash`'s main functionality.
Before proceeding, make sure [you have installed PaSh](../install/)

## Introduction

PaSh is a system for parallelizing POSIX shell scripts.
It has been shown to achieve order-of-magnitude performance improvements on shell scripts.

> _N.b.: PaSh is still under heavy development._

#### Example Script

Consider the following spell-checking script, applied to two large markdown files `f1.md` and `f2.md` (line 1):

```sh
# spell-checking.sh
cat f1.md f2.md | 
  tr A-Z a-z |
  tr -cs A-Za-z '\n' |
  sort |
  uniq | 
  comm -13 dict.txt - > out
cat out | wc -l | sed 's/$/ mispelled words!/'
```
The first `cat` streams two markdown files into a pipeline that converts characters in the stream into lower case, removes punctuation, sorts the stream in alphabetical order, removes duplicate words, and filters out words from a dictionary file (lines 1--7).
A second pipeline (line 7) counts the resulting lines to report the number of misspelled words to the user.

> If you're new to shell scripting, try to run each part of the pipeline separately and observe the output.
> For example, run `cat f1.md f2.md | tr A-Z a-z` in your terminal to witness the lower-case conversion.

Visually, the script can be thought as executing sequentially as follows:

<img src="https://docs.google.com/drawings/d/e/2PACX-1vQv-Krzb9hxWCbbQC9Zg5knm6SySJrayh3mdZXG3Z4Y6hC4kgQj4PWqYmxNAR-LyKN5Fu9lWHJV0J0F/pub?w=517&amp;h=55">

The first pipeline (left; parts omitted) processes `f1.md` and `f2.md` _sequentially_ through all pipeline stages and writes to `out`.
After it executes to completion, the second pipeline starts its _sequential_ execution.

#### Parallelizing Scripts with PaSh

PaSh transforms and executes each pipeline in a data-parallel fashion.
Visually, the parallel script would look like this for 2x-parallelism (i.e., assuming that the computer on which we execute the script has at least two CPUs and that PaSh is invoked with `-w` value of `2`).

<img src="https://docs.google.com/drawings/d/e/2PACX-1vR8AL-gkL7CvqRJYiyX8z20_WcJ68l9JUEinJgI-_jKussb6q33Qlc1saaXx7Cf2pYp8-qjKhYXGu5e/pub?w=517&amp;h=55">

Given a script, PaSh converts it to a dataflow graph, performs a series of semantics-preserving program transformations that expose parallelism, and then converts the dataflow graph back into a POSIX script.
The new parallel script has POSIX constructs added to explicitly guide parallelism, coupled with PaSh-provided Unix-aware runtime primitives for addressing performance- and correctness-related issues.

## Running Scripts

All scripts in this guide assume that `$PASH_TOP` is set to the top directory of the PaSh codebase (e.g., `/opt/pash` on docker)

**To run scripts in this section of the tutorial, make sure you are in the `intro` directory of the `evaluation`:**
```sh
cd $PASH_TOP/evaluation/intro
```

> In the following examples, you can avoid including `$PASH_TOP` before `pa.sh` by adding `PASH_TOP` in your `PATH`, which amounts to adding an `export PATH=$PATH:$PASH_TOP` in your shell configuration file.

#### Intro: Hello, Parallel World!

The simplest script to try out `pash` is `hello-world.sh`, which applies an expensive regular expression over the system's dictionary file.

To run `hello-world.sh` sequentially, you would call it using `bash`:
```sh
time bash ./hello-world.sh
```

To run it in parallel with PaSh:
```sh
time $PASH_TOP/pa.sh $PASH_TOP/evaluation/intro/hello-world.sh
```

At this point, you might be interested in running `pa.sh --help` to get a first sense of the available options.
Of particular interest is `--with` or `-w`, which specifies the degree of parallelism sought by PaSh (_e.g.,_ `-w 2`).

#### A More Interesting Script: Demo Spell

We will use `demo-spell.sh` -- a pipeline based [on the original Unix spell program](https://dl.acm.org/doi/10.1145/3532.315102) by Johnson -- to confirm that the infrastructure works as expected. We need to setup the appropriate input files for this script to execute:
```sh
./input/setup.sh
```
After inputs are configured, let's take a quick look at `demo-spell.sh`:
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
On our evaluation infrastructure, the 2x-parallel script takes about 28s.

You can check that the results are correct by:
```sh
diff spell.out pash-spell.out
```

Assuming you have more than 8 CPUs, you could also execute it with 8x-parallelism using:
```sh
time $PASH_TOP/pa.sh -w 8 -d 1 --log_file pash.log demo-spell.sh > pash-spell.out
``` 
On our evaluation infrastructure, the 8x-parallel script takes about 14s.

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

## What Next?

This concludes the first PaSh tutorial.
This section includes pointers for further exploration, depending on your needs.

#### The PaSh Repo

PaSh consist of three main components and a few additional "auxiliary" files and directories. 
The three main components are:

* [annotations](../../annotations/): DSL characterizing commands, parallelizability study, and associated annotations. More specifically, (i) a lightweight annotation language allows command developers to express key parallelizability properties about their commands; (ii) an accompanying parallelizability study of POSIX and GNU commands. guides the annotation language and optimized aggregator library 

* [compiler](../../compiler): Shell-dataflow translations and associated parallelization transformations. Given a script, the PaSh compiler converts it to a dataflow graph, performs a series of semantics-preserving program transformations that expose parallelism, and then converts the dataflow graph back into a POSIX script. 

* [runtime](../../runtime): Runtime components such as `eager`, `split`, and associated combiners. Apart from POSIX constructs added to guide parallelism explicitly, PaSh provides Unix-aware runtime primitives for addressing performance- and correctness-related issues.

These three components implement the contributions presented [in the EuroSys paper](https://arxiv.org/pdf/2007.09436.pdf).
They are expected to be usable with minimal effort, through a few different installation means presented below.

The auxiliary directories are:

* [docs](../../docs): Design documents, tutorials, installation instructions, etc.
* [evaluation](../../evaluation): Shell pipelines and script used for the evaluation of `pash`.

#### PaSh Concepts in Depth

Academic [papers](../README.md#academic-papers--events) associated with PaSh offer substantially deeper overviews of the concepts underpinning several PaSh components.

#### Useful Links

Chat:

* [Discord Server](https://discord.com/channels/947328962739187753/) ([Invite](http://join.binpa.sh/))

Mailing Lists: 

* [Discussion](https://groups.google.com/g/pash-dev): Join this mailing list for discussing all things `pash`
* [Commits](https://groups.google.com/g/pash-commits): Join this mailing list for commit notifications

Development/contributions:

* Contribution guide: [docs/contrib](../../docs/contributing/contrib.md)
* Continuous Integration Server: [http://ci.binpa.sh/](http://ci.binpa.sh/)

