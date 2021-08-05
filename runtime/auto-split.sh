#!/usr/bin/env bash

input="$1"
shift
outputs="$@"
n_outputs="$#"

# Set a default DISH_TOP in this directory if it doesn't exist
PASH_TOP=${PASH_TOP:-$(git rev-parse --show-toplevel)}

if type lsb_release >/dev/null 2>&1 ; then
   distro=$(lsb_release -i -s)
elif [ -e /etc/os-release ] ; then
   distro=$(awk -F= '$1 == "ID" {print $2}' /etc/os-release)
fi

# convert to lowercase
distro=$(printf '%s\n' "$distro" | LC_ALL=C tr '[:upper:]' '[:lower:]')
# now do different things depending on distro
case "$distro" in
   freebsd*)  
    temp="$(TMPDIR=/tmp mktemp -u pash_XXXXXXXXXX)"
    ;;
    *)
    temp="$(mktemp --tmpdir -u pash_XXXXXXXXXX)"
    ;;
esac


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
