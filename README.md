## PaSh: Light-touch Data-Parallel Shell Processing

Mailing lists: [Commits](mailto:pash-commits@googlegroups.com) | [Discussion](mailto:pash-discuss@googlegroups.com)

PaSh is  a system for  parallelizing POSIX shell  scripts. Key elements include:

* [compiler](./compiler): Shell-to-Dataflow translations and associated parallelization transformations.
* [annotations](./annotations/): DSL characterizing commands, parallelizability study, and associated annotations.
* [evaluation](./evaluation): shell pipelines and example [scripts](./evaluation/scripts) used for the evaluation.
* [runtime](./runtime): PaSh's runtime components, including `eager`, `split`, and assocaited combiners.
* [docs](./docs): Design documents, tutorials, installation instructions, etc.
* [papers](./papers): Academic papers related to PaSh ([EuroSys 2021](https://arxiv.org/abs/2007.09436)).

## Installation

If you are on a fresh clone, first run the following to download the `libdash` submodule:
```sh
git submodule init
git submodule update
```

Then run:
```sh
./install.sh -p
```

Don't forget to export the library path in the end :)
```sh
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/"
```

## Tests
To execute the current tests (which are all the one-liners) simply run the following:

```sh
cd compiler
./test_evaluation_scripts.sh
```

## Running PaSh

To parallelize script `./evaluation/hello-world.sh` with parallelization width of `2`, run:

```sh
./pa.sh -w 2 ./evaluation/hello-world.sh
``` 

