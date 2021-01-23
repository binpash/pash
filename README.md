## PaSh: Light-touch Data-Parallel Shell Processing
> PaSh is  a system for  parallelizing POSIX shell  scripts.

Quick Jump: [Repo Structure](#repo-structure) | [Running PaSh](#running-pash) | [Installation](#installation) | [Testing](#testing) | [Community & More](#community--more)

## Repo Structure

This repo hosts the core `pash` development. The structure is as follows:

* [annotations](./annotations/): DSL characterizing commands, parallelizability study, and associated annotations.
* [compiler](./compiler): Shell-Dataflow translations and associated parallelization transformations.
* [docs](./docs): Design documents, tutorials, installation instructions, etc.
* [evaluation](./evaluation): Shell pipelines and example [scripts](./evaluation/scripts) used for the evaluation.
* [runtime](./runtime): Runtime component — e.g., `eager`, `split`, and assocaited combiners.
* [scripts](./scripts): Scripts related to continuous integration, deployment, and testing.

## Running PaSh

To parallelize, say, `./evaluation/hello-world.sh` with parallelization width of `2`, from the top-level directory of the repository run:

```sh
./pa.sh ./evaluation/hello-world.sh
``` 

Run `./pa.sh --help` to get more information about the available commands.

## Installation

**Docker:** The easiest way to play with `pash` today is using Docker:

```sh
curl img.pash.ndr.md | docker load; docker run --name pash-playground -it pash/18.04
```

PaSh can be found in the container's `/pash` directory, so run `cd pash; git pull` to fetch the latest updates; more information in the [pash-on-docker guide](./docs/contrib.md#pash-on-docker-a-pocket-guide).

**Linux:** Alternatively, if you're on an Ubuntu, run:

```sh
curl up.pash.ndr.md | bash
```

## Tests

To execute the current tests, one-liner shell scripts, simply run:

```sh
cd compiler
./test_evaluation_scripts.sh
```

## Community & More

Mailing Lists: [Commits](https://groups.google.com/g/pash-commits) | [Discussion](https://groups.google.com/g/pash-discuss)

Continuous Integration Server: [http://pash.ndr.md/](http://pash.ndr.md/)

Contribution guide: [docs/contrib](docs/contrib.md)
