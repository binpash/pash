## PaSh: Light-touch Data-Parallel Shell Processing
> _A system for parallelizing POSIX shell scripts._

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
Read a longer tutorial, see [docs/tutorial](docs/tutorial.md).

## Installation

**Docker:** The easiest way to play with `pash` today is using Docker:

```sh
curl img.pash.ndr.md | docker load; docker run --name pash-playground -it pash/18.04
```

PaSh can be found in the container's `/pash` directory, so run `cd pash; git pull` to fetch the latest updates; more information in the [pash-on-docker guide](./docs/contrib.md#pash-on-docker-a-pocket-guide).

Alternatively, you can built the Docker container from scratch by running

```sh
git clone git@github.com:andromeda/pash.git
cd pash/scripts
docker build -t "pash/18.04" .
# Then you launch the container as above
docker run --name pash-playground -it pash/18.04
```

This should build a fresh docker image running `pash`.  If you wish to
run the tests, then you will need to extend the image.

```sh
cd pash/scripts
docker build -t pash-test --build-arg FROM_IMAGE=pash/18.04 - <Dockerfile.testing
docker run --name pash-with-tests -it pash-test
```

**Linux:** Alternatively, if you're on an Ubuntu 18.04, run:

```sh
curl up.pash.ndr.md | bash
```

This runs the install script in [scritps/install.sh](scritps/install.sh).
We have only tested this script on Ubuntu 18.04 (like the one used for the docker container and on AWS).

## Tests

To execute the current tests, one-liner shell scripts, simply run:

```sh
cd compiler
./test_evaluation_scripts.sh
```

## Community & More

Mailing Lists: 
* [Discussion](https://groups.google.com/g/pash-discuss): Join this mailing list for discussing all things `pash`
* [Commits](https://groups.google.com/g/pash-commits): Join this mailing list for commit notifications

Development/contributions:
* Contribution guide: [docs/contrib](docs/contrib.md)
* Continuous Integration Server: [ci.pash.ndr.md](http://ci.pash.ndr.md)
