#!/bin/bash

# trap ctrl-c and call ctrl_c()
trap cleanup INT

export DISH_TOP=${DISH_TOP:-${BASH_SOURCE%/*}}
export PASH_TOP=${PASH_TOP:-${DISH_TOP}/pash/}
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/"
# point to the local downloaded folders
export PYTHONPATH=${PASH_TOP}/python_pkgs/
export PASH_TIMESTAMP="$(date +"%y-%m-%d-%T")"

# Add hdfs support if hdfs command exist
if command -v "hdfs" &> /dev/null
then
    source "$PASH_TOP/compiler/dspash/hdfs_utils.sh"
fi

source "$PASH_TOP/compiler/orchestrator_runtime/pash_init_setup.sh" "$@" --distributed_exec

export PASH_TMP_PREFIX="$(mktemp -d /tmp/pash_XXXXXXX)/"

cleanup() {
        # kill "$FILEREADER_PID" "$DISCOVERY_PID"
        wait "$DISCOVERY_PID" 2>/dev/null
        rm -rf "$PASH_TMP_PREFIX"
}

# "$DISH_TOP/runtime/dspash/file_reader/filereader_server" &
# FILEREADER_PID=$!
"$DISH_TOP/runtime/bin/discovery_server" &
DISCOVERY_PID=$!
# python3 "$DISH_TOP/pash/compiler/dspash/worker.py" "$@" -d 1