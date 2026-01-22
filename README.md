## PaSh: Light-touch Data-Parallel Shell Processing

> _A system for parallelizing POSIX shell scripts._
> _Hosted by the [Linux Foundation](https://linuxfoundation.org/press-release/linux-foundation-to-host-the-pash-project-accelerating-shell-scripting-with-automated-parallelization-for-industrial-use-cases/)._

[![Tests](https://github.com/binpash/pash/actions/workflows/tests.yaml/badge.svg?branch=main&event=push)](https://github.com/binpash/pash/actions/workflows/tests.yaml?query=branch%3Amain++)
[![Build](https://github.com/binpash/pash/actions/workflows/build.yaml/badge.svg?branch=main&event=push)](https://github.com/binpash/pash/actions/workflows/build.yaml?query=branch%3Amain++)

Quick Jump: [Running PaSh](#running-pash) | [Installation](#installation) | [Testing](#testing) | [Repo Structure](#repo-structure) | [Community & More](#community--more) | [Citing](#citing)


## Installation via pip (Recommended)
```sh
pip install binpash-pash
pash -c "echo hi"
```

## Running PaSh

To parallelize, say, `./evaluation/intro/hello-world.sh` with parallelization degree of 2× run:

```sh
pash -w 2 ./evaluation/intro/hello-world.sh
```

If the script contains bash specific syntax, add the beta `--bash` flag to enable support.

Run `pash --help` to get more information about the available commands.
Jump to [docs/tutorial](docs/tutorial/) for a longer tutorial.


## Local testing (development)

To install and run PaSh for local development:

```sh
pip install -e .
pash --help
pash -c "echo hello | cat"
./scripts/run_tests.sh
```

For more details, manual installation, or other platforms see [installation instructions](./docs/install).

## Running with a local annotations library

To run with a local version of the library, please refer to the documentation [local annotations setup and usage](docs/local-annotations-library-documentation.md) 

## Repo Structure

This repo hosts the core `pash` development. The structure is as follows:

* [src/pash/](./src/pash/): Main Python package (installed via pip)
  * [preprocessor](./src/pash/preprocessor): Parses shell scripts, expands variables, and identifies dataflow regions for compilation.
  * [compiler](./src/pash/compiler): Translates shell dataflow regions to IRs and applies parallelization transformations.
  * [jit_runtime](./src/pash/jit_runtime): Just-in-time runtime that executes compiled regions and manages shell state.
  * [runtime](./src/pash/runtime): Runtime components — e.g., `eager`, `split`, and associated combiners.
* [docs](./docs): Design documents, tutorials, installation instructions, etc.
* [evaluation](./evaluation): Shell pipelines and example [scripts](./evaluation/other/more-scripts) used for the evaluation.
* [scripts](./scripts): Scripts related to continuous integration, deployment, and testing.

## Community & More

Chat:
* [Discord Server](ttps://discord.com/channels/947328962739187753/) ([Invite](https://discord.gg/6vS9TB97be))

Mailing Lists:
* [pash-devs](https://groups.google.com/g/pash-devs): Join this mailing list for discussing all things `pash`
* [pash-commits](https://groups.google.com/g/pash-commits): Join this mailing list for commit notifications

Development/contributions:
* Contribution guide: [docs/contributing](docs/contributing/contrib.md)
* Continuous Integration Server: [ci.binpa.sh](http://ci.binpa.sh)

## Citing

If you used PaSh, consider citing the following paper:
```
@inproceedings{pash2021eurosys,
author = {Vasilakis, Nikos and Kallas, Konstantinos and Mamouras, Konstantinos and Benetopoulos, Achilles and Cvetkovi\'{c}, Lazar},
title = {PaSh: Light-Touch Data-Parallel Shell Processing},
year = {2021},
isbn = {9781450383349},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
url = {https://doi.org/10.1145/3447786.3456228},
doi = {10.1145/3447786.3456228},
pages = {49–66},
numpages = {18},
keywords = {POSIX, Unix, pipelines, automatic parallelization, source-to-source compiler, shell},
location = {Online Event, United Kingdom},
series = {EuroSys '21}
}
```
