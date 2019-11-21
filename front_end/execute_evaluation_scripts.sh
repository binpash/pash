#!/bin/bash


./execute_evaluation_script.sh "minimal_sort_2" 
./execute_evaluation_script.sh "minimal_sort_4" 
./execute_evaluation_script.sh "minimal_sort_10" 
./execute_evaluation_script.sh "minimal_sort_20" 
./execute_evaluation_script.sh "minimal_sort_50" 
./execute_evaluation_script.sh "minimal_sort_100"

./execute_evaluation_script.sh "minimal_grep_1"
./execute_evaluation_script.sh "minimal_grep_2"
./execute_evaluation_script.sh "minimal_grep_4"
./execute_evaluation_script.sh "minimal_grep_10"
./execute_evaluation_script.sh "minimal_grep_20"
./execute_evaluation_script.sh "minimal_grep_50"
./execute_evaluation_script.sh "minimal_grep_100"
./execute_evaluation_script.sh "minimal_grep_200" 

./execute_evaluation_script.sh "spell_1" 

./execute_evaluation_script.sh "topn_2" 

./execute_evaluation_script.sh "wf_2" 
