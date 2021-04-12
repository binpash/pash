## PaSh: Light-touch Data-Parallel Shell Processing
> _A system for parallelizing POSIX shell scripts._

Quick Jump: [Running PaSh](#running-pash) | [Installation](#installation) | [Testing](#testing) | [Repo Structure](#repo-structure) | [Community & More](#community--more)

## Running PaSh

To parallelize, say, `./evaluation/intro/hello-world.sh` with parallelization width of `2`, from the top-level directory of the repository run:

```sh
./pa.sh ./evaluation/intro/hello-world.sh
``` 

Run `./pa.sh --help` to get more information about the available commands.
Read a longer tutorial, see [docs/tutorial](docs/tutorial.md).

## Installation

If on Ubuntu, Fedora, or Arch, you can simply run `curl up.pash.ndr.md | sh`. Alternatively, you can clone the repo and run `./scripts/up.sh`

You can also built the Docker container from scratch by running
```sh
git clone git@github.com:andromeda/pash.git
cd pash/scripts
docker build -t "pash/18.04" .
# Then launch continer:
docker run --name pash-play -it pash/18.04
```

More installation instructions in [the tutorial](./docs/tutorial.md#installation). Remember to `export` `PASH_TOP` in your startup scripts!

## Testing

To execute the current tests before committing or pushing code, simply run:

```sh
cd compiler
./test_evaluation_scripts.sh
```

More tests are available in [evaluation](./evaluation).

## Repo Structure

This repo hosts the core `pash` development. The structure is as follows:

* [annotations](./annotations/): DSL characterizing commands, parallelizability study, and associated annotations.
* [compiler](./compiler): Shell-Dataflow translations and associated parallelization transformations.
* [docs](./docs): Design documents, tutorials, installation instructions, etc.
* [evaluation](./evaluation): Shell pipelines and example [scripts](./evaluation/scripts) used for the evaluation.
* [runtime](./runtime): Runtime component — e.g., `eager`, `split`, and assocaited combiners.
* [scripts](./scripts): Scripts related to continuous integration, deployment, and testing.


## Community & More

Mailing Lists: 
* [Discussion](https://groups.google.com/g/pash-discuss): Join this mailing list for discussing all things `pash`
* [Commits](https://groups.google.com/g/pash-commits): Join this mailing list for commit notifications

Development/contributions:
* Contribution guide: [docs/contrib](docs/contrib.md)
* Continuous Integration Server: [ci.pash.ndr.md](http://ci.pash.ndr.md)
