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
To simply run PaSh on a script `script.sh` with parallelization width `2` make sure you are in the `compiler` directory and run:
```sh
python3.8 $PASH_TOP/compiler/pash.py --split_fan_out 2 script.sh
``` 

## Broader Impact

Ideas about broader impact can be found in this [document](https://docs.google.com/document/d/1XUAXr-Wt44Z2LLIN4OtK6FAlk-KOCHAs-_tWbKoJQGI/edit?usp=sharing).

