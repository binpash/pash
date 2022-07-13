## TODOs before merging to `future`

- working on all tests
- fix tests from compiler/test_evaluation_scripts.sh
  - alt_bigrams.sh
  - minimal_grep_stdin.sh: this one has different outputs when run with the script; however, if I add a cat <input> to it, it does not produce different inputs 
  - tr_test.sh: for that one, there is a line which reduces the whole data to one line and then the new model tries to RR-parallelize it still and manages to mess up things. An easy fix for that could be to basically disable parallelization after a tr which squeezes all new lines since there are no sequences of data to parallelize anyway for the moment. Long-term, we could allow parallelization but with a adj_line_merge aggregator. Outside the testing script, the outputs seem the same but somehow it still shows different outputs and I do not understand why.
  - shortest_scripts.sh: here I only needed to modify the script slightly: (1) option arguments for cut with whitespace as the parser cannot deal with them otherwise currently but we might want to change this in the future, (2) head -n 15 instead of head -15 which might be a bit harder to support and I did not really see how the man-page supports this actually when skimming but I might have missed that. 
- clean up utils for annotations
- graphviz
- Changing PaSh flags (making the default be priority r-split and then consecutive chunks), so remove the r_split flag and make defaults be the ones from the OSDI paper
- Fixing annotation library installation to a specific commit
- Remove code which got obsolete due to the changes