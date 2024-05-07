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

source "$PASH_TOP/compiler/pash_init_setup.sh" "$@" --distributed_exec

# export PASH_TMP_PREFIX="$(mktemp -d /tmp/pash_XXXXXXX)/"
export PASH_TMP_PREFIX="/tmp/pash_placeholder)"

function cleanup() {
        kill "$FILEREADER_PID" "$DISCOVERY_PID"
        wait "$FILEREADER_PID" "$DISCOVERY_PID" 2>/dev/null
        # rm -rf "$PASH_TMP_PREFIX"
}

"$DISH_TOP/runtime/bin/filereader_server" &
FILEREADER_PID=$!
"$DISH_TOP/runtime/bin/discovery_server" &
DISCOVERY_PID=$!

# No need to restart worker.py if we are resurrecting
# because worker.py is not killed in the first place!
if [[ "$@" != *"resurrect"* ]]; then
    python3 "$DISH_TOP/pash/compiler/dspash/worker.py" "$@"
fi
