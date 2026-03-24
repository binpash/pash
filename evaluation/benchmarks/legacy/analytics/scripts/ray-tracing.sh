#!/bin/bash
# source: posh benchmarks

in_dir=${1:-$IN}
out_dir=${2:-$OUT}
mkdir -p "$out_dir"

cat "$in_dir/1.INFO" | grep "\[RAY\]" | head -n1 | cut -c 7- > "$out_dir/rays.csv"
cat "$in_dir"/*.INFO | grep "\[RAY\]" | grep -v pathID | cut -c 7- >> "$out_dir/rays.csv"
cat "$out_dir/rays.csv" | sed -n '/^590432,/p' > "$out_dir/rt.log"


IN=analytics/inputs/ray_tracing_small/1.INFO

cat "$IN" | grep "\[RAY\]"  | cut -c 7- > "$out_dir/rays.csv"
cat "$out_dir/rays.csv" | sed -n '/^590432,/p' > "$out_dir/rt.log"



# cat "$in_dir/1.INFO" | grep "\[RAY\]" | head -n1 | cut -c 7- > "$out_dir/rays.csv"
# cat "$in_dir"/*.INFO | grep "\[RAY\]" | grep -v pathID | cut -c 7- >> "$out_dir/rays.csv"
cat "$out_dir/rays.csv" | sed -n '/^590432,/p' > "$out_dir/rt.log"