#!/bin/bash

source "$PASH_TOP/compiler/orchestrator_runtime/pash_orch_lib.sh"

export PASH_SPEC_NODE_DIRECTORY="${PASH_TMP_PREFIX}/speculative/partial_order/"

pash_spec_communicate_scheduler()
{
    local message=$1
    pash_communicate_unix_socket "PaSh-Spec-scheduler" "${PASH_SPEC_SCHEDULER_SOCKET}" "${message}"
}

pash_spec_communicate_scheduler_just_send()
{
    pash_spec_communicate_scheduler "$1"
}

pash_spec_wait_until_scheduler_listening()
{
    pash_wait_until_unix_socket_listening "PaSh-Spec-scheduler" "${PASH_SPEC_SCHEDULER_SOCKET}"
}


start_server()
{
    python3 -S "$PASH_SPEC_TOP/parallel-orch/scheduler_server.py" "$@" &
    export daemon_pid=$!
    ## Wait until daemon has established connection
    pash_spec_wait_until_scheduler_listening
}

cleanup_server()
{
    local daemon_pid=$1
    ## Only wait for daemon if it lives (it might be dead, rip)
    if ps -p "$daemon_pid" > /dev/null
    then
        ## Send and receive from daemon
        msg="Done"
        daemon_response=$(pash_spec_communicate_scheduler "$msg")
        wait 2> /dev/null 1>&2 
    fi
}

export -f pash_spec_communicate_scheduler
export -f pash_spec_communicate_scheduler_just_send
export -f pash_spec_wait_until_scheduler_listening
export -f start_server
export -f cleanup_server

