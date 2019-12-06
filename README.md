# Dish: Automatically-Distributed Shell Scripting

The repository contains:

* [compiler](./compiler): Dish's data-flow model and associated optimization passes.
* [distributability](./distributability/): Distributability analysis used to inform the DSL characterizing commands.
* [evaluation](./evaluation): shell pipelines and example [scripts](./evaluation/scripts) used for the evaluation.
* [libdash](./libdash): Git submodule pointing to the `libdash` POSIX-compliant parser.
* [parser](./parser): The parser uses the `libdash` POSIX-compliant parser and `atdgen` to convert scripts to JSON-encoded ASTs.
* [pldi](./pldi): The paper submitted to PLDI 2020.
* [runtime](./runtime): A JVM-based runtime for Dish.

