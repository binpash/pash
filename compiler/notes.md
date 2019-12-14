## Supported Scripts

The following files just contain POSIX (no process substitution) pipes
and commands

- compile.sh.json" - OK - xargs
- grep.sh.json" # OK - single input file
- grep2.sh.json" # OK - 2 input files
- minimal.sh.json" # ISSUE - redirection of a file - have to be handled in the translation
- minimal2.sh.json" # OK - 2 input files
- minimal3.sh.json" # OK - 1 input file
- minimal4.sh.json" # OK - 1 input file
- minimal5.sh.json" # OK - 1 input file
- shortest-scripts.sh.json" # PROBLEM - unexpanded *
- spell.sh.json" # PROBLEM - sort and others
- topn.sh.json" # PROBLEM - sort and others
- wc.sh.json" # OK - single input file
- wf.sh.json" # PROBLEM - sort and others
- ngrams.sh.json" # PROBLEM - sort, tee, and others
- diff.sh.json" # PROBLEM - Function definition is not handled. Can this be solved?

The following contain And operators together with pipes and commands

- old_ngrams.sh.json" # ISSUE - redirection of a file - have to be handled in the translation

The following is interesting, since it contains command substitution
(which is parsed as backticks in the libdash parser). The backticks
seem to mean that whatever is in the backticks will execute first, and
its output will become a string in place of the backticks

- page-count.sh.json"

Unidentified:

- json_filename="../scripts/json/maximal.sh.json"
