#!/usr/bin/env bash

export PASH_TOP=${PASH_TOP:-${BASH_SOURCE%/*}}
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/"
# point to the local downloaded folders
export PYTHONPATH=${PASH_TOP}/python_pkgs/

if [ "$#" -eq 1 ] && [ "$1" = "--init" ]; then
  $PASH_TOP/compiler/superoptimize.sh
  exit
fi

if ! command -v python3 &> /dev/null
then
    echo "Python >=3 could not be found"
    exit
fi

## Create a temporary directory where PaSh can use for temporary files and logs
export PASH_TMP_PREFIX="$(mktemp -d /tmp/pash_XXXXXXX)/"

## Create the input and output fifo that the runtime will use for communication
export RUNTIME_IN_FIFO="${PASH_TMP_PREFIX}/runtime_in_fifo"
export RUNTIME_OUT_FIFO="${PASH_TMP_PREFIX}/runtime_out_fifo"
## TODO: Get rid of these two commands if possible
rm -f "$RUNTIME_IN_FIFO" "$RUNTIME_OUT_FIFO"
mkfifo "$RUNTIME_IN_FIFO" "$RUNTIME_OUT_FIFO"
export DAEMON_SOCKET="${PASH_TMP_PREFIX}/daemon_socket"

python3 -S "$PASH_TOP/compiler/pash_runtime_daemon.py" $@ &
daemon_pid=$!

## Initialize all things necessary for pash to execute (logging/functions/etc)
source "$PASH_TOP/compiler/pash_init_setup.sh" "$@"

## Wait until daemon has established connection
##
## TODO: Can we get rid of the `sleep` in this wait?
pash_wait_until_daemon_listening

PASH_FROM_SH="pa.sh" python3 -S $PASH_TOP/compiler/pash.py "$@"
pash_exit_code=$?

## Only wait for daemon if it lives (it might be dead, rip)
if ps -p $daemon_pid > /dev/null 
then
  ## Send and receive from daemon
  msg="Done"
  daemon_response=$(pash_communicate_daemon "$msg")
  wait 2> /dev/null 1>&2 
fi

## Don't delete the temporary directory if we are debugging
if [ "$PASH_DEBUG_LEVEL" -eq 0 ]; then
  rm -rf "${PASH_TMP_PREFIX}"
fi

(exit $pash_exit_code)
