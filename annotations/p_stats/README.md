The study has evolved over time

Originally, there are five main groups:

* _stateless_: these are the simplest to parallelize. They fall into four classes depending on the input chunk---some are parallelizable at the level of individual characters (e.g., `tr`), some at the level of lines (e.g., `grep`), some at the level of paragraphs (e.g., fmt), and some at the level of files

* _pure_: these are somewhat more difficult to parallelize, as they might require to see the end of the input (or include some line metadata); they are still pure (a la reducers).

* _DFS_: these interact with the (distributed version of the) file-system. As the file-system is a central part of the Unix design and philosophy, many of these commands have Unix-specific semantics that are not helpful in pipelines or distributed workloads---e.g., operations related to file ownership. However, they are not difficult to emulate atop a conventional distributed storage.

* _EVs_: these commands affect environment variables; interestingly, the vast majority of these commands only _read_ environment variables---the common case in scripts.
So even if we had something like distributed transactions, we would still be able to get away with mostly read-accesses that do not need to block (note that there are multiple transactional protocols, because these are reversible).

* _Side-effectful_: These have non-reversible side-effects that affect the outside system. They are non-distributable, and would most certainly require a pesimistic transaction management control.

There are still a few commands whose semantics are not entirely understood, marked as `?` (but for which `pash` assumes the most conservative class possible).


Currently, the four major classes:

* _Stateless Commands_: This class contains commands that operate on individual line elements of their input, without maintaining state across invocations.
These are commands that can be expressed as a purely functional `map` or `filter`—e.g., `grep` filters out individual lines and `basename` removes a path prefix from a string.
They may produce multiple elements—e.g., `tr` may insert NL tokens—but always return empty output for empty input.
Workloads that use only stateless commands are trivial to parallelize: they do not require any synchronization to maintain correctness, nor caution about where to split inputs.

* _Parallelizable Pure Commands_: This class contains commands that respect functional purity—i.e., same outputs for same inputs—but maintain internal state across their entire pass.
The details of this state and its propagation during element processing affect their parallelizability characteristics.
Some commands are easy to parallelize, because they maintain trivial state and are commutative—e.g., `wc` simply maintains a counter.
Other commands, such as `sort`, maintain more complex invariants that have to be taken into account when merging partial results.

* _Non-parallelizable Pure Commands_: This class contains commands that, while purely functional, cannot be parallelized.
This is because their internal state depends on prior state in the same pass in non-trivial ways.
For example, hashing commands such as `sha1sum` maintain complex state that has to be updated sequentially.
If parallelized on a single input, each stage would need to wait on the results of all previous stages, foregoing any parallelism benefits.

* _Side-effectful Commands_: This class contains commands that have side-effects across the system—for example, updating environment variables, interacting with the
filesystem, and accessing the network.
Such commands are not parallelizable without finer-grained concurrency control mechanisms that can detect side-effects across the system.
This is the largest class.
