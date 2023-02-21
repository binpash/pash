#!/bin/bash

# trap ctrl-c and call ctrl_c()
trap cleanup INT

export PASH_TOP=${PASH_TOP:-${BASH_SOURCE%/*}}
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/"
# point to the local downloaded folders
export PYTHONPATH=${PASH_TOP}/python_pkgs/
export PASH_TIMESTAMP="$(date +"%y-%m-%d-%T")"

# add hdfs directory if hdfs command exist
if command -v "hdfs" &> /dev/null
then
    datanode_dir=$(hdfs getconf -confKey dfs.datanode.data.dir) 
    export HDFS_DATANODE_DIR=${datanode_dir#"file://"} # removes file:// prefix
fi

source "$PASH_TOP/compiler/orchestrator_runtime/pash_init_setup.sh" "$@" --distributed_exec

export PASH_TMP_PREFIX="$(mktemp -d /tmp/pash_XXXXXXX)/"

cleanup() {
        kill "$FILEREADER_PID" "$DISCOVERY_PID"
        wait "$FILEREADER_PID" "$DISCOVERY_PID" 2>/dev/null
        rm -rf "$PASH_TMP_PREFIX"
}

"$PASH_TOP/runtime/dspash/file_reader/filereader_server" &
FILEREADER_PID=$!
"$PASH_TOP/runtime/dspash/file_reader/discovery_server" &
DISCOVERY_PID=$!
python3 "$PASH_TOP/compiler/dspash/worker.py" "$@"
