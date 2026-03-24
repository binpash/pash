#!/bin/bash

BENCHMARK_DIR=$(cd "$(dirname "$0")" && pwd)
input_dir=$BENCHMARK_DIR/inputs
URL=https://atlas.cs.brown.edu/data
mkdir -p "$input_dir"

size=full
for arg in "$@"; do
    case "$arg" in
    --small) size=small ;;
    --min) size=min ;;
    esac
done

if [ -f "$input_dir/routeviews.mrt" ]; then
    echo "routeviews.mrt already exists, skipping download."
else
    wget --no-check-certificate "${URL}/log-analysis/routeviews.mrt" -O "$input_dir/routeviews.mrt"
fi

if [[ "$size" == "min" ]]; then
    if [[ ! -d "$input_dir/ray_tracing_$size" ]]; then
        wget --no-check-certificate $URL/log-analysis/ray_tracing_$size.tar.gz -O "$input_dir/ray_tracing_$size.tar.gz"
        tar -xzf "$input_dir/ray_tracing_$size.tar.gz" -C "$input_dir" --no-same-owner
        rm "$input_dir/ray_tracing_$size.tar.gz"
    fi
    if [[ ! -d "$input_dir/pcaps_$size" ]]; then
        mkdir -p "$input_dir/pcaps_$size"
        cp "${BENCHMARK_DIR}/min_inputs/pcaps/"* "$input_dir/pcaps_$size"
    fi
    if [[ ! -d "$input_dir/nginx-logs_$size" ]]; then
        mkdir -p "$input_dir/nginx-logs_$size"
        cp "${BENCHMARK_DIR}/min_inputs/nginx-logs/"* "$input_dir/nginx-logs_$size"
    fi
    if [[ ! -d "$input_dir/port_scan_$size" ]]; then
        wget --no-check-certificate $URL/log-analysis/port_scan_$size.tar.gz -O "$input_dir/port_scan_$size.tar.gz"
        tar -xzf "$input_dir/port_scan_$size.tar.gz" -C "$input_dir" --no-same-owner
        rm "$input_dir/port_scan_$size.tar.gz"
    fi
    exit 0
elif [[ "$size" == "small" ]]; then
    if [[ ! -d "$input_dir/ray_tracing_$size" ]]; then
        wget --no-check-certificate $URL/log-analysis/ray_tracing_$size.tar.gz -O "$input_dir/ray_tracing_$size.tar.gz"
        tar -xzf "$input_dir/ray_tracing_$size.tar.gz" -C "$input_dir" --no-same-owner
        rm "$input_dir/ray_tracing_$size.tar.gz"
    fi
    if [[ ! -d "$input_dir/pcaps_$size" ]]; then
        wget --no-check-certificate $URL/pcaps.zip -O "$input_dir/pcaps_$size.zip"
        unzip "$input_dir/pcaps_$size.zip" -d "$input_dir"
        mv "$input_dir/pcaps" "$input_dir/pcaps_$size"
        rm "$input_dir/pcaps_$size.zip"
    fi
    if [[ ! -d "$input_dir/nginx-logs_$size" ]]; then
        zip_dst="$input_dir/nginx.zip"
        wget --no-check-certificate $URL/nginx.zip -O "$zip_dst"
        unzip "$zip_dst" -d "$input_dir"
        mv "$input_dir/nginx-logs" "$input_dir/nginx-logs_$size"
        rm "$zip_dst"
    fi
    if [[ ! -d "$input_dir/port_scan_$size" ]]; then
        wget --no-check-certificate $URL/log-analysis/port_scan_$size.tar.gz -O "$input_dir/port_scan_$size.tar.gz"
        tar -xzf "$input_dir/port_scan_$size.tar.gz" -C "$input_dir" --no-same-owner
        rm "$input_dir/port_scan_$size.tar.gz"
    fi
    exit 0
else
    if [[ ! -d "$input_dir/ray_tracing_$size" ]]; then
        wget --no-check-certificate $URL/log-analysis/ray_tracing_$size.tar.gz -O "$input_dir/ray_tracing_$size.tar.gz"
        tar -xzf "$input_dir/ray_tracing_$size.tar.gz" -C "$input_dir" --no-same-owner
        rm "$input_dir/ray_tracing_$size.tar.gz"
    fi
    if [[ ! -d "$input_dir/pcaps_$size" ]]; then
        zip_dst="$input_dir/pcaps.zip"
        wget --no-check-certificate "$URL"/pcaps.zip -O "$zip_dst"
        unzip "$zip_dst" -d "$input_dir"
        mv "$input_dir/pcaps/" "$input_dir/pcaps_$size/"
        rm "$zip_dst"

        wget --no-check-certificate "$URL"/pcaps_large.zip -O "$zip_dst"
        unzip "$zip_dst" -d "$input_dir"/pcaps_$size
        rm "$zip_dst"
    fi
    if [[ ! -d "$input_dir/nginx-logs_$size" ]]; then
        zip_dst="$input_dir/nginx.zip"
        wget --no-check-certificate $URL/nginx.zip -O "$zip_dst"
        unzip "$zip_dst" -d "$input_dir"
        mv "$input_dir/nginx-logs" "$input_dir/nginx-logs_$size"
        rm "$zip_dst"

        zip_dst="$input_dir/nginx_large.zip"
        wget --no-check-certificate "$URL"/log-analysis/web-server-access-logs.zip -O "$zip_dst"
        unzip "$zip_dst" -d "$input_dir"
        mv "$input_dir/access.log" "$input_dir/nginx-logs_$size/access.log"
        rm "$zip_dst"
    fi
    if [[ ! -d "$input_dir/port_scan_$size" ]]; then
        wget --no-check-certificate $URL/log-analysis/port_scan_$size.tar.gz -O "$input_dir/port_scan_$size.tar.gz"
        tar -xzf "$input_dir/port_scan_$size.tar.gz" -C "$input_dir" --no-same-owner
        rm "$input_dir/port_scan_$size.tar.gz"
    fi
fi
