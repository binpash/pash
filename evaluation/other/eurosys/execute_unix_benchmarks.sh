#!/bin/bash

## Necessary to set PASH_TOP
export PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel --show-superproject-working-tree)}

## This sets up to what extent we run the evaluation.
## There are 2 levels:
## 1. Small inputs (1GB) | --width 4
## 2. Big inputs (10GB) | --width 16 (EuroSys evaluation)
evaluation_level=1

while getopts 'slh' opt; do
    case $opt in
        s) evaluation_level=1 ;;
        l) evaluation_level=2 ;;
        h) echo "There are two possible execution levels:"
           echo "option -s: Small inputs (1GB) | --width 4"
           echo "option -l: Big inputs (10GB) | --width 16 (EuroSys evaluation)"
           exit 0 ;;
        *) echo 'Error in command line parsing' >&2
           exit 1
    esac
done
shift "$(( OPTIND - 1 ))"


unixfun_dir="$PASH_TOP/evaluation/unixfun/"
unixfun_intermediary="${unixfun_dir}/intermediary/"
intermediary_dir="$PASH_TOP/evaluation/intermediary/"
results_subdir_prefix="unixfun"

if [ "$evaluation_level" -eq 1 ]; then
    echo "Executing Unix50 scripts with 1GB inputs and --width 4"
    maximum_input_size="$((1024 * 1024 * 1024))" # 1 GB
    n_in=4
elif [ "$evaluation_level" -eq 2 ]; then
    echo "Executing Unix50 scripts with 10GB inputs and --width 16"
    maximum_input_size="$((10 * 1024 * 1024 * 1024))" # 10 GB
    n_in=16
else
    echo "Unrecognizable execution level: $evaluation_level"
    exit 1
fi

results_subdir="${results_subdir_prefix}_${n_in}_${maximum_input_size}"

rm -r $unixfun_intermediary
mkdir -p $unixfun_intermediary
mkdir -p $intermediary_dir
mkdir -p "$PASH_TOP/evaluation/results/${results_subdir}/"

## Make inputs larger and generate scripts and their envs
python3 generate_unixfun_scripts.py $unixfun_dir $unixfun_intermediary $maximum_input_size

for unixfun_pipeline in $(ls ${unixfun_intermediary} | grep -v "_env" | cut -f 1 -d '.' | sort); do
    echo $unixfun_pipeline

    echo "Generating input and intermediary scripts... be patient..."
    python3 "$PASH_TOP/evaluation/generate_microbenchmark_intermediary_scripts.py" \
            $unixfun_intermediary $unixfun_pipeline $n_in $intermediary_dir

    echo "Executing script with bash and pash..."
    "$PASH_TOP/evaluation/execute_compile_evaluation_script.sh" -s -a "${unixfun_pipeline}" "${n_in}" "${results_subdir}" > /dev/null 2>&1
    rm -f /tmp/eager*
done
