#!/bin/bash


./execute_compile_evaluation_script.sh "minimal_sort_2" $distr_output_dir 1 0
./execute_compile_evaluation_script.sh "minimal_sort_4" $distr_output_dir 1 0
./execute_compile_evaluation_script.sh "minimal_sort_10" $distr_output_dir 1 0
./execute_compile_evaluation_script.sh "minimal_sort_20" $distr_output_dir 1 0
./execute_compile_evaluation_script.sh "minimal_sort_50" $distr_output_dir 1 0
./execute_compile_evaluation_script.sh "minimal_sort_100" $distr_output_dir 1 0

./execute_compile_evaluation_script.sh "minimal_grep_2" $distr_output_dir 1 0
./execute_compile_evaluation_script.sh "minimal_grep_4" $distr_output_dir 1 0
./execute_compile_evaluation_script.sh "minimal_grep_10" $distr_output_dir 1 0
./execute_compile_evaluation_script.sh "minimal_grep_20" $distr_output_dir 1 0
./execute_compile_evaluation_script.sh "minimal_grep_50" $distr_output_dir 1 0
./execute_compile_evaluation_script.sh "minimal_grep_100" $distr_output_dir 1 0
./execute_compile_evaluation_script.sh "minimal_grep_200" $distr_output_dir 1 0

# ./execute_compile_evaluation_script.sh "spell_1" 

## This has to be tuned (fan-out and batch size)
./execute_compile_evaluation_script.sh "topn_2" $distr_output_dir 1 0
./execute_compile_evaluation_script.sh "topn_4" $distr_output_dir 1 0
./execute_compile_evaluation_script.sh "topn_10" $distr_output_dir 1 0
./execute_compile_evaluation_script.sh "topn_20" $distr_output_dir 1 0
./execute_compile_evaluation_script.sh "topn_50" $distr_output_dir 1 0
./execute_compile_evaluation_script.sh "topn_100" $distr_output_dir 1 0

./execute_compile_evaluation_script.sh "wf_2" $distr_output_dir 1 0
./execute_compile_evaluation_script.sh "wf_4" $distr_output_dir 1 0
./execute_compile_evaluation_script.sh "wf_10" $distr_output_dir 1 0
./execute_compile_evaluation_script.sh "wf_20" $distr_output_dir 1 0
./execute_compile_evaluation_script.sh "wf_50" $distr_output_dir 1 0
./execute_compile_evaluation_script.sh "wf_100" $distr_output_dir 1 0
