## PaSh: Light-touch Data-Parallel Shell Processing

> _A system for parallelizing POSIX shell scripts._
> _Hosted by the [Linux Foundation](https://linuxfoundation.org/press-release/linux-foundation-to-host-the-pash-project-accelerating-shell-scripting-with-automated-parallelization-for-industrial-use-cases/)._


| Service      | Main | Develop |
| :---         |    :----:   |          :----: |
| Tests      | [![Tests](https://github.com/binpash/pash/actions/workflows/tests.yaml/badge.svg?branch=main&event=push)](https://github.com/binpash/pash/actions/workflows/tests.yaml?query=branch%3Amain++)     | [![Tests](https://github.com/binpash/pash/actions/workflows/tests.yaml/badge.svg?branch=future&event=push)](https://github.com/binpash/pash/actions/workflows/tests.yaml?query=branch%3Afuture++)   |
| Build   | [![Build](https://github.com/binpash/pash/actions/workflows/build.yaml/badge.svg?branch=main&event=push)](https://github.com/binpash/pash/actions/workflows/build.yaml?query=branch%3Amain++)        | [![Build](https://github.com/binpash/pash/actions/workflows/build.yaml/badge.svg?branch=future&event=push)](https://github.com/binpash/pash/actions/workflows/build.yaml?query=branch%3Afuture++)      |
| Pages     | [![DeployPublish](https://github.com/binpash/web/actions/workflows/pages.yaml/badge.svg)](https://github.com/binpash/web/actions/workflows/pages.yaml) | [![DeployPublish](https://github.com/binpash/web/actions/workflows/pages.yaml/badge.svg)](https://github.com/binpash/web/actions/workflows/pages.yaml) |


Quick Jump: [Running PaSh](#running-pash) | [Installation](#installation) | [Testing](#testing) | [Repo Structure](#repo-structure) | [Community & More](#community--more) | [Citing](#citing)

## Running PaSh

To parallelize, say, `./evaluation/intro/hello-world.sh` with parallelization degree of 2× run:

```sh
./pa.sh ./evaluation/intro/hello-world.sh
```

Run `./pa.sh --help` to get more information about the available commands.
Jump to [docs/tutorial](docs/tutorial/) for a longer tutorial.

## Installation

On Ubuntu, Fedora, Debian, or Arch, run `curl up.pash.ndr.md | sh` to set up PaSh.

For more details, manual installation, or other platforms see [installation instructions](./docs/install).

## Testing

To execute the current tests before committing and pushing code, simply run:

```sh
./scripts/run_tests.sh
```

## Repo Structure

This repo hosts the core `pash` development. The structure is as follows:

* [annotations](./annotations/): DSL characterizing commands, parallelizability study, and associated annotations.
* [compiler](./compiler): Shell-dataflow translations and associated parallelization transformations.
* [docs](./docs): Design documents, tutorials, installation instructions, etc.
* [evaluation](./evaluation): Shell pipelines and example [scripts](./evaluation/other/more-scripts) used for the evaluation.
* [runtime](./runtime): Runtime component — e.g., `eager`, `split`, and associated combiners.
* [scripts](./scripts): Scripts related to continuous integration, deployment, and testing.

## Community & More

Chat:
* [Discord Server](ttps://discord.com/channels/947328962739187753/) ([Invite](http://join.binpa.sh/))

Mailing Lists:
* [pash-dev](https://groups.google.com/g/pash-dev): Join this mailing list for discussing all things `pash`
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

Caruca notes:
sudo docker build -t pash:caruca -f scripts/docker/debian/Dockerfile .

sudo docker run -it -d -v "$(pwd)":/opt/pash --name pash_caruca pash:caruca bash
sudo docker exec -it pash_caruca bash

Run all benchmarks: evaluation/osdi22-eval/run_all.sh
