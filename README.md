# Dish: Automatically-Distributed Shell Scripting

The repository contains:

* [compiler](./compiler): Dish's data-flow model and associated optimization passes.
* [distributability](./distributability/): Distributability analysis used to inform the DSL characterizing commands.
* [evaluation](./evaluation): shell pipelines and example [scripts](./evaluation/scripts) used for the evaluation.
* [parser](./parser): The parser uses the `libdash` POSIX-compliant parser and `atdgen` to convert scripts to JSON-encoded ASTs.
* [pldi](./pldi): The paper submitted to PLDI 2020.
* [runtime](./runtime): A JVM-based runtime for Dish.
* [scripts](./scripts): Auxiliary scripts
* [tests](./tests): Contains a set of tests to check that everything works

# High level Overview

Figure containing high level overview of the Dish workflow:
https://docs.google.com/drawings/d/15f3JKa-xF7_yyumX2YGyOWOFmZO_KBog1QC1PQr_LJw/edit

# Some TODOs

Minor things left out by the PLDI submission:

* Collect results for two additional _macro_-benchmarks: (i) a bio-informatics pipeline, and (ii) web-indexing pipeline.
* Report single-node performance results for all macro-benchmarks

Additions for OSDI:

* Implement, and possibley extend the DSL---and ensure it works with existing commands

* POSIX Compliance: scripts have to be able to run with any POSIX script, even if they do not speed up much of it. That is, there shouldn't be a POSIX script that fails, only ones that do not see any performance improvements.

* Environment variables: in the distributed execution, add the ability to synchronize environment variables (and other environment context) across nodes. This requires defining a bit more precisely what makes sense to distribute (e.g., PIDs and machine names cannot be distributed, and Dish should be careful with `/proc/fs`).

* Distributed file-system (and, more generally, networked distribution): as much of Unix is about files, incorporating an existing Unix-like distributed file-system is key to shell script distribution. With a DFS, `/x/y` would resolve to the filesystem of a specific node---but it should be location-aware so that Dish sends commands to where files (streams) reside.

* Planner: augment and extend the functionality of the planner with the ability to optimize the distributed script and map it efficiently on the underlying infrastructure---using information about the workload, infrastructure, etc. The planner should also add and leverage metadata (lines, batches etc.) for key planning decisions.

* JIT (partial evaluation): introduce the ability to call the planner at runtime, when the latest information possible becomes available. Example of information that is not visible includes values of variables (that can be streams), files, shell expansions etc.

* Distributed runtime environment: currently Dish rewrites pipelines to extract parallelism---that is, the end result of our toolchain is still a shell pipeline. A more sophisticated runtime environment could provide significant benefits by manipulating metadata for both accelerating performance (e.g., wrapping commands with packing/unpacking) and fault tolerance (by knowing which of calls fail). It is important for the runtime to be optimized significantly so as to not have any issues.
