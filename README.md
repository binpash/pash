## PaSh: Light-touch Data-Parallel Shell Processing

Mailing lists: [Commits](mailto:pash-commits@googlegroups.com) ([join](https://groups.google.com/g/pash-commits)) | [Discussion](mailto:pash-discuss@googlegroups.com) ([join](https://groups.google.com/g/pash-discuss))

PaSh is  a system for  parallelizing POSIX shell  scripts. Key elements include:

* [compiler](./compiler): Shell-to-Dataflow translations and associated parallelization transformations.
* [annotations](./annotations/): DSL characterizing commands, parallelizability study, and associated annotations.
* [evaluation](./evaluation): shell pipelines and example [scripts](./evaluation/scripts) used for the evaluation.
* [runtime](./runtime): PaSh's runtime components, including `eager`, `split`, and assocaited combiners.
* [docs](./docs): Design documents, tutorials, installation instructions, etc.
* [papers](./papers): Academic papers related to PaSh ([EuroSys 2021](https://arxiv.org/abs/2007.09436)).

## Running PaSh

To parallelize, say, `./evaluation/hello-world.sh` with parallelization width of `2`, from the top-level directory of the repository run:

```sh
./pa.sh -w 2 ./evaluation/hello-world.sh
``` 

## Installation

**Docker:** The easiest way to play with `pash` today is using Docker:

```sh
curl -s img.pash.ndr.md | docker load; docker run -i pash-latest
```

After you're in the image, run `cd pash; git pull` to get the latest updates.

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

