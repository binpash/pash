#!/bin/bash

unix50_dir="../evaluation/unix50/"
unix50_intermediary="${unix50_dir}/intermediary/"
intermediary_dir="../evaluation/intermediary/"

mkdir -p $unix50_intermediary
mkdir -p $intermediary_dir

## Make inputs larger and generate scripts and their envs
input_size_increase=10
python3 generate_unix50_scripts.py $unix50_dir $unix50_intermediary $input_size_increase

## TODO: Generate intermediaries

## TODO: Execute them
