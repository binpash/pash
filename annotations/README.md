# Parallelizability Study & Annotation Language
Quick Jump: [Parallelizability](#main-parallelizability-classes) | [study](#parallelizability-study-of-commands-in-gnu--posix) | [example 1](#a-simple-example-chmod) | [example 1](#another-example-cut) | [howto](#how-to-annotate-a-command) | [issues](#Issues)

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
The second class, `pure`, contains commands that respect functional purity -- _i.e.,_ same outputs for same inputs -- but maintain internal state across their entire pass.
The details of this state and its propagation during element processing affect their parallelizability characteristics.
Some commands are easy to parallelize, because they maintain trivial state and are commutative -- _e.g.,_ `wc` simply maintains a counter.
Other commands, such as `sort`, maintain more complex invariants that have to be taken into account when merging partial results.

* _Non-parallelizable Pure Commands:_
The third class, `non-parallelizable`, contains commands that, while purely functional, cannot be parallelized within a single data stream.
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

The first step to annotating a command is to identify its default class: `stateless`, `pure`, `non-parallelizable`, and `side-effectful`. How does the command behave without any inputs?
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

## Issues

* [pr.json line 18] (https://github.com/binpash/pash/tree/main/annotations/pr.json)
* [cat.json line 3] (https://github.com/binpash/pash/tree/main/annotations/cat.json)

