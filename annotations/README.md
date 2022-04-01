# Parallelizability Study & Annotation Language
Quick Jump: [Parallelizability](#main-parallelizability-classes) | [study](#parallelizability-study-of-commands-in-gnu--posix) | [example 1](#a-simple-example-chmod) | [example 2](#another-example-cut) | [howto](#how-to-annotate-a-command) | [issues](#Issues)

PaSh includes 
  (i) a parallelizability study of commands in POSIX and GNU Coreutils, and 
  (ii) an annotation language for describing the parallelizability properties of individual commands.
The parallelizability study informed the design of the annotation language, which was in turn used to capture the key parallelizability characteristics in many of these commands.

> _N.b.: We welcome contributions to the study and annotatations for common commands._

## Main Parallelizability Classes

PaSh introduces four major parallelizability classes:

* _Stateless Commands:_
The first class, `stateless`, contains commands that operate on individual line elements of their input, without maintaining state across invocations.
These are commands that can be expressed as a purely functional `map` or `filter` -- _e.g.,_ `grep` filters out individual lines and `basename` removes a path prefix from a string.
They may produce multiple elements -- _e.g.,_ `tr` may insert `NL` tokens -- but always return empty output for empty input.
Workloads that use only stateless commands are trivial to parallelize:
  they do not require any synchronization to maintain correctness, nor caution about where to split inputs.

* _Parallelizable Pure Commands:_
The second class, `parallelizable_pure`, contains commands that respect functional purity -- _i.e.,_ same outputs for same inputs -- but maintain internal state across their entire pass.
The details of this state and its propagation during element processing affect their parallelizability characteristics.
Some commands are easy to parallelize, because they maintain trivial state and are commutative -- _e.g.,_ `wc` simply maintains a counter.
Other commands, such as `sort`, maintain more complex invariants that have to be taken into account when merging partial results.

* _Non-parallelizable Pure Commands:_
The third class, `pure`, contains commands that, while purely functional, cannot be parallelized within a single data stream.
This is because their internal state depends on prior state in non-trivial ways over the same pass. % should we say something about state machines?
For example, hashing commands such as `sha1sum` maintain complex state that has to be updated sequentially.
If parallelized on a single input, each stage would need to wait on the results of all previous stages, foregoing any parallelism benefits.

* _Side-effectful Commands:_
The last class, `side-effectful`, contains commands that have side-effects across the system -- for example, updating environment variables, interacting with the filesystem, and accessing the network.
Such commands are not parallelizable without finer-grained concurrency control mechanisms that can detect side-effects across the system.

## Parallelizability Study of Commands in GNU & POSIX

The parallelizability study of commands in GNU and POSIX is comprised of two parts: a coarse-grained parallelizability study, and a set of annotations for commands.

The main results of the parallelizability study are summarized in the [PaSh EuroSys'21 paper (Sec. 3.1 and Tab. 1)](https://arxiv.org/pdf/2007.09436.pdf).
To see the results of the  parallelizability study, run [./p_stats](./p_stats).

Annotations for about 60 popular commands are stored in this directory encoded as JSON files (about 14 lines per annotation on average, for a total of 846 lines of annotations).
Annotations can be thought of as defining a bidirectional correspondence between a command and a node in the dataflow graph---the abstraction used by the PaSh compiler.
Since command behaviors (and correspondence) can change based on their arguments, annotations contain a sequence of predicates.
Each predicate is accompanied by information that instantiates the correspondence between a command and a dataflow node.

## A Simple Example: `chmod`

As a first example, below we present the annotations for `chmod`.

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

## Another Example: `cut`

As another example, below we present the annotations for `cut`.

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
The first predicate indicates that `cut` is "pure" -- _i.e._, not parallelizable but representable as a dataflow node -- if the value accompanying the `-d` option is `\n` or if it was used with the `-z` flag.
In both of these cases, newlines do not represent data item boundaries, but are rather used internally by the command, making it unsafe to parallelize by splitting on line boundaries.
In all other cases (see the "default" case) the command is stateless.
Inputs are always assigned to the non-option arguments and the output is always stdout.
The option "stdin-hyphen" indicates that a non-option argument that is just a dash `-` represents the stdin, and the option “empty-args-stdin” indicates that if non-option arguments are empty, then the command reads from its stdin.
The list identified by "short-long" contains a correspondence of short and long argument names for this command.

## How to Annotate a Command

The first step to annotating a command is to identify its default class: `stateless`, `parallelizable_pure`, `pure`, and `side-effectful`. How does the command behave without any inputs?
The next step is to identify the set of inputs and their order.

This process then has to be repeated for every set of arguments, which have to be expressed as first-order-logic predicates (see examples above).
This can be (and is currently) achieved in an incremental fashion:
  a few flags at a time.

For more details, here is an early version of the annotation language:

```
  <option> ::= `-' <string>
  <category> ::= `stateless' | `pure' | ...
  <maybe-int> ::= ε | <int>
  <arg> ::= `args[' <int> `]'
  <args> ::= <arg>
           | `args[' <maybe-int> `:' <maybe-int> `]'
  <input> ::= `stdin' 
            |  <args>
  <inputs> ::= <input>
             | <input> `,' <inputs>
  <output> ::= `stdout' 
             |  <arg>
  <outputs> ::= <output>
              | <output> `,' <outputs>
  <option-pred> ::= <option>
                  | `value' <option> = <string>
                  | `not' <option-pred>
                  | <option-pred> `or' <option-pred>
                  | <option-pred> `and' <option-pred>
  <assignment> ::= `(' <category>, `[' <inputs> `]' `,' `[' <output> `]' `)'
  <predicate> ::= <option-pred> `=>' <assignment>
  <pred-list> ::= `|' <predicate> <pred-list>
                | `|' `otherwise' `=>' <assignment>
  <command> ::= <name> `\{' <pred-list> `\}'
  <command-list> ::= <command>
                   | <command> <command-list>
```

[//]: # (TODO: 1. update language spec; 2. put all annotations in a directory)

## Mini-tutorial: Adding Custom Aggregators

For this tutorial, let's assume you want to parallelize [a simple `ann-agg.sh` script](https://github.com/binpash/pash/blob/main/evaluation/tests/ann-agg.sh).

Let's also assume there are no annotations or aggregators for the commands `test_one` and `test_two`.
Note that normally these two commands would be annotated as `stateless`, as their aggregator is simply the con`cat`enation function;
  however, we will now annotate them as `parallelizable_pure` and provide "custom" aggregation commands that simply concatenate their input streams.

*Step 1: Implement aggregators and their annotations*:

An aggregator is usually either binary or _n_-ary:
  it takes as input two or _n_ file names (or paths) and outputs results to the standard out.
An aggregator  may also take additional flags---for example, flags that configure its operation or flags that were provided to the original command.

We will implement `test_one`'s aggregator as [a shell script](https://github.com/binpash/pash/blob/main/runtime/agg/opt/concat.sh) that internally uses the Unix `cat` command to concatenate any number of input streams.

We will implement `test_two`'s aggregator as [a Python script](https://github.com/binpash/pash/blob/main/runtime/agg/py/cat.py) that concatenates any number of inputs streams.

For PaSh to be able to hook these aggregators correctly, _i.e._, so that it can instantiate them as command invocations, we also need to add their annotations in [annotations/custom_aggregators](https://github.com/binpash/pash/tree/main/annotations/custom_aggregators).
Below are the two annotation files named [`annotations/custom_aggregators/cat.py.json`](./custom_aggregators/cat.py.json) and [`annotations/custom_aggregators/concat.json`](./custom_aggregators/concat.json). (FIXME: relative path? **Until this is fixed, prefix aggregator names with `pagg-` to avoid name clashes!**)
The most important information in these files is (i) the aggregation command's `name`, and (ii) its treatment of inputs (both taking `["args[:]"]`), and outputs (both outputing to `["stdout"]`).

*Step 2: Point commands to their custom aggregators*:
Add two new annotation files in `$PASH_TOP/annotations` with names `test_one.json` and  `test_two.json`, so that they point to the right aggregator commands.
Apart from providing the correct command `name`, the two key properties are the `class` (which should be `parallelizable_pure`) and the `rel_path` (which should point to the aggregator programs we just implemented---ideally, relative to `$PASH_TOP`).

Here is the annotation for [`test_one.json`](./test_one.json), where the aggregator points to `runtime/agg/opt/concat.sh`. 
Note that path is relative with respect to `$PASH_TOP` and therefore refers to `$PASH_TOP/runtime/agg/opt/concat.sh`:

Here is the annotation for [`test_two.json`](./test_two.json), pointing to `runtime/agg/py/cat.py` (i.e., implying `$PASH_TOP/runtime/agg/py/cat.py`).
The annotations also specifies that the aggregator should be called with the `-a` flag, in addition to any other flags provided to the original command.

**More complex aggregators**:
Suppose we want to parallelize a new script called [ann-agg-2.sh](https://github.com/binpash/pash/blob/main/evaluation/tests/ann-agg.sh).
This script contains two new commands `test_uniq_1` and `test_uniq_2`. 
Their annotations are in files [annotations/test_uniq_1](./test_uniq_1.json) and [annotations/test_uniq_2.json](./test_uniq_2.json).

## Issues

* [pr.json line 18] (https://github.com/binpash/pash/tree/main/annotations/pr.json)
* [cat.json line 3] (https://github.com/binpash/pash/tree/main/annotations/cat.json)

