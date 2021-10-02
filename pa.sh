#!/usr/bin/env bash

export PASH_TOP=${PASH_TOP:-${BASH_SOURCE%/*}}
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/"

if [ "$#" -eq 1 ] && [ "$1" = "--init" ]; then
  $PASH_TOP/compiler/superoptimize.sh
  exit
fi

if ! command -v python3 &> /dev/null
then
    echo "Python >=3 could not be found"
    exit
fi

## Get distro
## TODO: Move that somewhere where it happens once (during installation)
if type lsb_release >/dev/null 2>&1 ; then
    distro=$(lsb_release -i -s)
elif [ -e /etc/os-release ] ; then
    distro=$(awk -F= '$1 == "ID" {print $2}' /etc/os-release)
fi

# convert to lowercase
export distro=$(printf '%s\n' "$distro" | LC_ALL=C tr '[:upper:]' '[:lower:]')

## Create a temporary directory where PaSh can use for temporary files and logs
export PASH_TMP_PREFIX="$(mktemp -d /tmp/pash_XXXXXXX)/"

## TODO: Add these fifos in the random directory
export RUNTIME_IN_FIFO=runtime_in_fifo
export RUNTIME_OUT_FIFO=runtime_out_fifo
rm -f "$RUNTIME_IN_FIFO" "$RUNTIME_OUT_FIFO"
mkfifo "$RUNTIME_IN_FIFO" "$RUNTIME_OUT_FIFO"

python3 "$PASH_TOP/compiler/pash_runtime_daemon.py" "$RUNTIME_IN_FIFO" "$RUNTIME_OUT_FIFO" $@ &
daemon_pid=$!

PASH_FROM_SH="pa.sh" python3 $PASH_TOP/compiler/pash.py "$@"

## TODO: Make sure to properly terminate the pash runtime "daemon"
kill -15 "$daemon_pid" 2> /dev/null 1>&2
wait 2> /dev/null 1>&2 
