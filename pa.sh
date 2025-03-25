#!/usr/bin/env bash

export PASH_TOP=${PASH_TOP:-${BASH_SOURCE%/*}}
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/"
# point to the local downloaded folders
export PYTHONPATH="${PASH_TOP}/python_pkgs/:${PYTHONPATH}"
## Register the signal handlers, we can add more signals here
trap kill_all SIGTERM SIGINT

## kill all the pending processes that are spawned by this shell
kill_all() {
    # kill all my subprocesses only
    kill -s SIGKILL 0
    # kill pash_daemon
    kill -s SIGKILL "$daemon_pid"
}
## Save the umask to first create some files and then revert it
old_umask=$(umask)

## Restore the umask to create files etc
umask u=rwx,g=rx,o=rx

if [ "$#" -eq 1 ] && [ "$1" = "--init" ]; then
  "$PASH_TOP"/compiler/superoptimize.sh
  exit
fi

if ! command -v python3 &> /dev/null
then
    echo "Python >=3 could not be found"
    exit
fi

## get bash version for pash
export PASH_BASH_VERSION="${BASH_VERSINFO[@]:0:3}"

## Create a temporary directory where PaSh can use for temporary files and logs
export PASH_TMP_PREFIX="$(mktemp -d /tmp/pash_XXXXXXX)/"

## Create a timestamp that PaSh can use for log directories 
##   (should not be used to create critical directories/files, only logs/monitors/etc,
##    all the cricial pash temp files should go in PASH_TMP_PREFIX)
export PASH_TIMESTAMP="$(date +"%y-%m-%d-%T")"

## Create the input and output fifo that the runtime will use for communication
export RUNTIME_IN_FIFO="${PASH_TMP_PREFIX}/runtime_in_fifo"
export RUNTIME_OUT_FIFO="${PASH_TMP_PREFIX}/runtime_out_fifo"
## TODO: Get rid of these two commands if possible
rm -f "$RUNTIME_IN_FIFO" "$RUNTIME_OUT_FIFO"
mkfifo "$RUNTIME_IN_FIFO" "$RUNTIME_OUT_FIFO"
export DAEMON_SOCKET="${PASH_TMP_PREFIX}/daemon_socket"
export DSPASH_SOCKET="${PASH_TMP_PREFIX}/dspash_socket"

## Initialize all things necessary for pash to execute (logging/functions/etc)
source "$PASH_TOP/compiler/orchestrator_runtime/pash_init_setup.sh" "$@"

## This starts a different server depending on the configuration
if [ "$show_version" -eq 0 ]; then
  ## Exports $daemon_pid
  start_server "$@"
fi

## Restore the umask before executing
umask "$old_umask"
PASH_FROM_SH="pa.sh" python3 -S "$PASH_TOP/compiler/pash.py" "$@"
pash_exit_code=$?
if [ "$show_version" -eq 0 ]; then
  cleanup_server "${daemon_pid}"
fi

## Don't delete the temporary directory if we are debugging
if [ "$PASH_DEBUG_LEVEL" -eq 0 ]; then
  rm -rf "${PASH_TMP_PREFIX}"
fi

(exit "$pash_exit_code")
