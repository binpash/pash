#!/bin/bash

distr_output_dir=/tmp/dish_output/

## Call the translate for all the scripts
./translate_script.sh "../evaluation/minimal_sort_2" $distr_output_dir 1 0
./translate_script.sh "../evaluation/minimal_sort_4" $distr_output_dir 1 0
./translate_script.sh "../evaluation/minimal_sort_10" $distr_output_dir 1 0
./translate_script.sh "../evaluation/minimal_sort_20" $distr_output_dir 1 0
./translate_script.sh "../evaluation/minimal_sort_50" $distr_output_dir 1 0
./translate_script.sh "../evaluation/minimal_sort_100" $distr_output_dir 1 0

./translate_script.sh "../evaluation/minimal_grep_1" $distr_output_dir 1 0
./translate_script.sh "../evaluation/minimal_grep_2" $distr_output_dir 1 0
./translate_script.sh "../evaluation/minimal_grep_4" $distr_output_dir 1 0
./translate_script.sh "../evaluation/minimal_grep_10" $distr_output_dir 1 0
./translate_script.sh "../evaluation/minimal_grep_20" $distr_output_dir 1 0
./translate_script.sh "../evaluation/minimal_grep_50" $distr_output_dir 1 0
./translate_script.sh "../evaluation/minimal_grep_100" $distr_output_dir 1 0
./translate_script.sh "../evaluation/minimal_grep_200" $distr_output_dir 1 0

./translate_script.sh "../evaluation/spell_1" $distr_output_dir

./translate_script.sh "../evaluation/topn_2" $distr_output_dir 2 200000
./translate_script.sh "../evaluation/topn_4" $distr_output_dir 4 200000
./translate_script.sh "../evaluation/topn_10" $distr_output_dir 10 200000
./translate_script.sh "../evaluation/topn_20" $distr_output_dir 20 200000
./translate_script.sh "../evaluation/topn_50" $distr_output_dir 50 200000
./translate_script.sh "../evaluation/topn_100" $distr_output_dir 100 200000


./translate_script.sh "../evaluation/wf_2" $distr_output_dir
