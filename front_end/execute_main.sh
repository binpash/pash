#!/bin/bash

# The following files just contain POSIX (no process substitution)
# pipes and commands
#
# json_filename="../scripts/json/compile.sh.json" # OK - xargs
# json_filename="../scripts/json/grep.sh.json" # OK - single input file
# json_filename="../scripts/json/grep2.sh.json" # OK - 2 input files
# json_filename="../scripts/json/minimal.sh.json" # ISSUE - redirection of a file - have to be handled in the translation
# json_filename="../scripts/json/minimal2.sh.json" # OK - 2 input files
# json_filename="../scripts/json/minimal3.sh.json" # OK - 1 input file
# json_filename="../scripts/json/minimal4.sh.json" # OK - 1 input file
json_filename="../scripts/json/minimal5.sh.json" # OK - 1 input file
# json_filename="../scripts/json/shortest-scripts.sh.json" # PROBLEM - unexpanded *
# json_filename="../scripts/json/spell.sh.json" # PROBLEM - sort and others
# json_filename="../scripts/json/topn.sh.json" # PROBLEM - sort and others
# json_filename="../scripts/json/wc.sh.json" # OK - single input file
# json_filename="../scripts/json/wf.sh.json" # PROBLEM - sort and others
# json_filename="../scripts/json/ngrams.sh.json" # PROBLEM - sort, tee, and others
# json_filename="../scripts/json/diff.sh.json" # PROBLEM - Function definition is not handled. Can this be solved?

# The following contain And operators together with pipes and commands
#
# json_filename="../scripts/json/old_ngrams.sh.json" # ISSUE - redirection of a file - have to be handled in the translation


# The following is interesting, since it contains command substitution
# (which is parsed as backticks in the Greenberg parser). The
# backticks seem to mean that whatever is in the backticks will
# execute first, and its output will become a string in place of the
# backticks
# json_filename="../scripts/json/page-count.sh.json"

# Unidentified
#
# json_filename="../scripts/json/maximal.sh.json"

python3 main.py $json_filename
