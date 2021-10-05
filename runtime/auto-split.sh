#!/usr/bin/env bash

input="$1"
shift
outputs="$@"
n_outputs="$#"

# Set a default DISH_TOP in this directory if it doesn't exist
PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}
# TODO change all auto-split calls to provide $distro as argument
if type lsb_release >/dev/null 2>&1 ; then
   distro=$(lsb_release -i -s)
elif [ -e /etc/os-release ] ; then
   distro=$(awk -F= '$1 == "ID" {print $2}' /etc/os-release)
fi

# convert to lowercase
distro=$(printf '%s\n' "$distro" | LC_ALL=C tr '[:upper:]' '[:lower:]')
# generate a temp file
temp="$(TMPDIR=/tmp mktemp -t pash_XXXXXXXXXX)"

cat "$input" > "$temp"
total_lines=$(wc -l $temp | cut -f 1 -d ' ')
batch_size=$( expr $total_lines / $n_outputs )
# echo "Input: $input"
# echo "Ouputs: $outputs"
# echo "Number of outputs: $n_outputs"
# echo "Total Lines: $total_lines"
# echo "Batch Size: $batch_size"

cleanup()
{
    kill -SIGPIPE $split_pid > /dev/null 2>&1
}
trap cleanup EXIT


# echo "$PASH_TOP/evaluation/tools/split $input $batch_size $outputs"
$PASH_TOP/runtime/split "$temp" "$batch_size" $outputs &
split_pid=$!
wait $split_pid
rm -f $temp
