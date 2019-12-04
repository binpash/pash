Currently, there are five main groups:

* _stateless_: these are the simplest to parallelize. They fall into four classes depending on the input chunk---some are parallelizable at the level of individual characters (e.g., `tr`), some at the level of lines (e.g., `grep`), some at the level of paragraphs (e.g., fmt), and some at the level of files

* _pure_: these are somewhat more difficult to parallelize, as they might require to see the end of the input (or include some line metadata); they are still pure (a la reducers).

* _DFS_: these interact with the (distributed version of the) file-system. As the file-system is a central part of the Unix design and philosophy, many of these commands have Unix-specific semantics that are not helpful in pipelines or distributed workloads---e.g., operations related to file ownership. However, they are not difficult to emulate atop a conventional distributed storage.

* _EVs_: these commands affect environment variables; interestingly, the vast majority of these commands only _read_ environment variables---the common case in scripts.
So even if we had something like distributed transactions, we would still be able to get away with mostly read-accesses that do not need to block (note that there are multiple transactional protocols, because these are reversible).

* _Side-effectful_: These have non-reversible side-effects that affect the outside system. They are non-distributable, and would most certainly require a pesimistic transaction management control.

There are still a few commands whose semantics are not entirely understood, marked as `?` (but which should fall under one of these categories).


Raflib
