#!/bin/bash

input="$1"
shift
outputs="$@"
n_outputs="$#"

# Set a default DISH_TOP in this directory if it doesn't exist
DISH_TOP=${DISH_TOP:-$(git rev-parse --show-toplevel)}

temp=$(mktemp -u)

cat "$input" > "$temp"
total_lines=$(wc -l $temp | cut -f 1 -d ' ')
batch_size=$( expr $total_lines / $n_outputs )
# echo "Input: $input"
# echo "Ouputs: $outputs"
# echo "Number of outputs: $n_outputs"
# echo "Total Lines: $total_lines"
# echo "Batch Size: $batch_size"
# echo "$DISH_TOP/evaluation/tools/split $input $batch_size $outputs"
$DISH_TOP/evaluation/tools/split "$temp" "$batch_size" $outputs
