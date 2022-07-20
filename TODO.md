## TODOs before merging to `future`

- fix tests from compiler/test_evaluation_scripts.sh:
  + bigrams
- clean up utils for annotations
- graphviz
- Changing PaSh flags (making the default be priority r-split and then consecutive chunks), so remove the r_split flag and make defaults be the ones from the OSDI paper
- Fixing annotation library installation to a specific commit
- Remove code which got obsolete due to the changes
- Room for optimization: basically disable parallelization after a tr which squeezes all new lines since there are no sequences of data to parallelize anyway for the moment. 
    Long-term, we could allow parallelization but with a adj_line_merge aggregator.
- Changes to scripts:
  + `spell-grep.sh`: added `export -f set_diff` as the command was not known when parallelizing?
  + `shortest_scripts.sh`: here I only needed to modify the script slightly: 
    (1) option arguments for `cut` with whitespace as the parser cannot deal with them otherwise currently but we might want to change this in the future, 
    (2) `head -n 15` instead of `head -15` which might be a bit harder to support. I did not really see how the man-page supports this actually when skimming but I might have missed that. 
- tr_test.sh: Outside the testing script, the outputs are the same but somehow it still shows different outputs. Checked this with Konstantinos and he will check the testing script later.
