#!/bin/bash

# pash_input_args=(1 2 "3 4" "%")
# echo "${pash_input_args[@]@Q}"
printf "\"%s\" " "${pash_input_args[@]}"

