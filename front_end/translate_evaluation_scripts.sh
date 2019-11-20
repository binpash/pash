#!/bin/bash

distr_output_dir=/tmp/dish_output/

## Call the translate for all the scripts
./translate_script.sh "../evaluation/minimal_sort_2" $distr_output_dir
./translate_script.sh "../evaluation/minimal_sort_4" $distr_output_dir
./translate_script.sh "../evaluation/minimal_sort_10" $distr_output_dir
./translate_script.sh "../evaluation/minimal_sort_20" $distr_output_dir
./translate_script.sh "../evaluation/minimal_sort_50" $distr_output_dir
./translate_script.sh "../evaluation/minimal_sort_100" $distr_output_dir
