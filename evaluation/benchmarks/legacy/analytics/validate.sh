#!/bin/bash

TOP=$(git rev-parse --show-toplevel)
eval_dir="${TOP}/evaluation/benchmarks/analytics"
hashes_dir="${eval_dir}/hashes"
outputs_dir="${eval_dir}/outputs"
mkdir -p "${outputs_dir}"

size="full"
generate=false
selected_scripts=""

while [ $# -gt 0 ]; do
    case "$1" in
        --generate)
            generate=true
            shift
            ;;
        --small)
            size="small"
            shift
            ;;
        --min)
            size="min"
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

cd "$outputs_dir" || exit

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

if $generate; then
    if should_run "pcaps"; then
        md5sum "pcaps_$size"/* > "$hashes_dir/pcaps_$size.md5sum"
    fi

    if should_run "nginx"; then
        md5sum "nginx_$size"/* > "$hashes_dir/nginx_$size.md5sum"
    fi

    if should_run "port-scan"; then
        md5sum port_scan_$size/as_popularity.csv > "$hashes_dir/port_scan_$size.md5sum"
    fi

    if should_run "ray-tracing"; then
        md5sum ray_tracing_$size/* > "$hashes_dir/ray_tracing_$size.md5sum"
    fi

    exit 0
fi

if should_run "pcaps"; then
    bench=pcaps_$size
    md5sum --check --quiet --status "$hashes_dir/$bench.md5sum"
    echo $bench $?
fi

if should_run "nginx"; then
    bench=nginx_$size
    md5sum --check --quiet --status "$hashes_dir/$bench.md5sum"
    echo $bench $?
fi

if should_run "port-scan"; then
    bench=port_scan_$size
    md5sum --check --quiet --status "$hashes_dir/$bench.md5sum"
    echo $bench $?
fi

if should_run "ray-tracing"; then
    bench=ray_tracing_$size
    md5sum --check --quiet --status "$hashes_dir/$bench.md5sum"
    echo $bench $?
fi
