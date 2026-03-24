#!/bin/bash

TOP=$(git rev-parse --show-toplevel)
eval_dir="${TOP}/evaluation/benchmarks/analytics"
input_dir="${eval_dir}/inputs"
scripts_dir="${eval_dir}/scripts"
outputs_dir="${eval_dir}/outputs"
mkdir -p "$outputs_dir"

export LC_ALL=C

size=full
selected_scripts=""

while [ $# -gt 0 ]; do
    case "$1" in
        --small)
            size=small
            shift
            ;;
        --min)
            size=min
            shift
            ;;
        -s|--scripts)
            shift
            while [ $# -gt 0 ] && [ "$(echo "$1" | cut -c1)" != "-" ]; do
                if [ -z "$selected_scripts" ]; then
                    selected_scripts="$1"
                else
                    selected_scripts="$selected_scripts $1"
                fi
                shift
            done
            ;;
        *)
            shift
            ;;
    esac
done

nginx_input=$input_dir/nginx-logs_$size
pcaps_input=$input_dir/pcaps_$size
port_scan_input=$input_dir/port_scan_$size/all_logs.jsonl
rt_inputs_dir=$input_dir/ray_tracing_$size
rt_outputs_dir=$outputs_dir/ray_tracing_$size

export BENCHMARK_CATEGORY="analytics"
KOALA_SHELL=${KOALA_SHELL:-bash}

should_run() {
    script_name=$1
    # If no scripts specified, run all
    if [ -z "$selected_scripts" ]; then
        return 0
    fi
    for selected in $selected_scripts; do
        if [ "$selected" = "$script_name" ]; then
            return 0
        fi
    done
    return 1
}

# nginx benchmark
if should_run "nginx"; then
    echo "nginx"
    BENCHMARK_INPUT_FILE="$(realpath "$nginx_input")"
    export BENCHMARK_INPUT_FILE
    BENCHMARK_SCRIPT="$(realpath "$scripts_dir/nginx.sh")"
    export BENCHMARK_SCRIPT

    $KOALA_SHELL $scripts_dir/nginx.sh $nginx_input $outputs_dir/nginx_$size
    echo $?
fi

# pcaps benchmark
if should_run "pcaps"; then
    echo "pcaps"
    BENCHMARK_INPUT_FILE="$(realpath "$pcaps_input")"
    export BENCHMARK_INPUT_FILE

    BENCHMARK_SCRIPT="$(realpath "$scripts_dir/pcaps.sh")"
    export BENCHMARK_SCRIPT

    $KOALA_SHELL $scripts_dir/pcaps.sh $pcaps_input $outputs_dir/pcaps_$size
    echo $?
fi

# port-scan benchmark
if should_run "port-scan"; then
    echo "port-scan"
    BENCHMARK_INPUT_FILE="$(realpath "$port_scan_input")"
    export BENCHMARK_INPUT_FILE

    BENCHMARK_SCRIPT="$(realpath "$scripts_dir/port-scan.sh")"
    export BENCHMARK_SCRIPT
    go_install_dir="${eval_dir}/go_install"

    export PATH=$PATH:/$go_install_dir/go/bin
    export GOPATH=$HOME/go
    export PATH=$PATH:$GOPATH/bin

    mkdir -p "$outputs_dir/port_scan_$size"
    touch "$outputs_dir/port_scan_$size/annotated.jsonl"
    touch "$outputs_dir/port_scan_$size/file1"
    touch "$outputs_dir/port_scan_$size/file2"
    touch "$outputs_dir/port_scan_$size/as_popularity.csv"

    $KOALA_SHELL $scripts_dir/port-scan.sh "$port_scan_input" "$eval_dir/inputs/routeviews.mrt" "$outputs_dir/port_scan_$size/annotated" "$outputs_dir/port_scan_$size/file1" "$outputs_dir/port_scan_$size/file2" "$outputs_dir/port_scan_$size/as_popularity"
    echo $?
fi

# ray-tracing benchmark
if should_run "ray-tracing"; then
    echo "ray-tracing"
    mkdir -p "$outputs_dir/ray_tracing_$size"
    BENCHMARK_SCRIPT="$(realpath "$scripts_dir/ray-tracing.sh")"
    export BENCHMARK_SCRIPT
    BENCHMARK_INPUT_FILE="$(realpath "$rt_inputs_dir")"
    export BENCHMARK_INPUT_FILE
    $KOALA_SHELL "$scripts_dir/ray-tracing.sh" "$rt_inputs_dir" "$rt_outputs_dir"
    echo $?
fi
