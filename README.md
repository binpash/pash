## PaSh: Light-touch Data-Parallel Shell Processing

Mailing lists: [Commits](https://groups.google.com/g/pash-commits) | [Discussion](https://groups.google.com/g/pash-discuss)

PaSh is  a system for  parallelizing POSIX shell  scripts. Key elements include:

* [annotations](./annotations/): DSL characterizing commands, parallelizability study, and associated annotations.
* [compiler](./compiler): Shell-Dataflow translations and associated parallelization transformations.
* [docs](./docs): Design documents, tutorials, installation instructions, etc.
* [evaluation](./evaluation): Shell pipelines and example [scripts](./evaluation/scripts) used for the evaluation.
* [papers](./papers): Academic papers related to PaSh ([EuroSys 2021](https://arxiv.org/abs/2007.09436)).
* [runtime](./runtime): Runtime component — e.g., `eager`, `split`, and assocaited combiners.
* [scripts](./scripts): Scripts related to continuous integration, deployment, and testing.

## Running PaSh

To parallelize, say, `./evaluation/hello-world.sh` with parallelization width of `2`, from the top-level directory of the repository run:

```sh
./pa.sh -w 2 ./evaluation/hello-world.sh
``` 

## Installation

**Docker:** The easiest way to play with `pash` today is using Docker:

```sh
curl img.pash.ndr.md | docker load; docker run --name pash-playground -it pash/18.04
```

PaSh can be found in the container's `/pash` directory, so run `cd pash; git pull` to fetch the latest updates; more information in the [pash-on-docker guide](./docs/contrib.md#pash-on-docker-a-pocket-guide).

**Linux:** Alternatively, if you're on an Ubuntu, run:

```sh
curl -s up.pash.ndr.md | bash
```

Appending `-- -a` additionally installs dependencies such as `opam`, `python3`, etc. but requires `sudo` (i.e., "root"). This option is great for AWS images.

## Tests

To execute the current tests, one-liner shell scripts, simply run:

```sh
cd compiler
./test_evaluation_scripts.sh
```

