# source the local pash config
[ -f ~/.pash_init ] && source ~/.pash_init
## File directory
export RUNTIME_DIR=$(dirname "${BASH_SOURCE[0]}")
export WRAPPER_LIB_DIR="$RUNTIME_DIR/../wrapper_library/"
## TODO: Is there a better way to do this?
export RUNTIME_LIBRARY_DIR="$RUNTIME_DIR/../../runtime/"
export PASH_REDIR="&2"
export PASH_DEBUG_LEVEL=0

## Check flags
export pash_output_time_flag=1
export pash_execute_flag=1
export pash_dry_run_compiler_flag=0
export pash_assert_compiler_success_flag=0
export pash_assert_all_regions_parallelizable_flag=0
export pash_checking_log_file=0
export pash_checking_debug_level=0
export pash_avoid_pash_runtime_completion_flag=0
export pash_profile_driven_flag=1
export pash_no_parallel_pipelines=0
export pash_daemon_communicates_through_unix_pipes_flag=0
export pash_speculative_flag=0
export show_version=0
export distributed_exec=0

for item in "$@"
do

    if [ "$pash_checking_log_file" -eq 1 ]; then
        export pash_checking_log_file=0
        export PASH_REDIR="$item"
    fi

    if [ "$pash_checking_debug_level" -eq 1 ]; then
        export pash_checking_debug_level=0
        export PASH_DEBUG_LEVEL=$item
    fi

    # We output time always 
    # if [ "--output_time" == "$item" ]; then
    #     pash_output_time_flag=1
    # fi
    if [ "--version" == "$item" ]; then
        export show_version=1
    fi
    if [ "--dry_run_compiler" == "$item" ]; then
        export pash_dry_run_compiler_flag=1
    fi

    if [ "--assert_compiler_success" == "$item" ]; then
        export pash_assert_compiler_success_flag=1
    fi

    if [ "--assert_all_regions_parallelizable" == "$item" ]; then
        export pash_assert_all_regions_parallelizable_flag=1
    fi

    if [ "--log_file" == "$item" ]; then
        pash_checking_log_file=1
    fi

    if [ "--avoid_pash_runtime_completion" == "$item" ]; then
        export pash_avoid_pash_runtime_completion_flag=1
    fi

    if [ "--profile_driven" == "$item" ]; then
        export pash_profile_driven_flag=1
    fi

    if [ "-d" == "$item" ] || [ "--debug" == "$item" ]; then
        pash_checking_debug_level=1
    fi

    if [ "--no_parallel_pipelines" == "$item" ]; then
        export pash_no_parallel_pipelines=1
    fi

    if [ "--daemon_communicates_through_unix_pipes" == "$item" ]; then
        export pash_daemon_communicates_through_unix_pipes_flag=1
    fi

    if [ "--speculative" == "$item" ]; then
        export pash_speculative_flag=1
    fi

    if [ "--distributed_exec" == "$item" ]; then
        export distributed_exec=1
    fi
done

## `pash_redir_output` and `pash_redir_all_output` are strictly for logging.
##
## They do not execute their arguments if there is no debugging.
if [ "$PASH_DEBUG_LEVEL" -eq 0 ]; then
    pash_redir_output()
    {
        :
    }

    pash_redir_all_output()
    {
        :
    }

    pash_redir_all_output_always_execute()
    {
        > /dev/null 2>&1 "$@"
    }

else
    if [ "$PASH_REDIR" == '&2' ]; then
        pash_redir_output()
        {
            >&2 "$@"
        }

        pash_redir_all_output()
        {
            >&2 "$@"
        }

        pash_redir_all_output_always_execute()
        {
            >&2 "$@"
        }
    else
        pash_redir_output()
        {
            >>"$PASH_REDIR" "$@"
        }

        pash_redir_all_output()
        {
            >>"$PASH_REDIR" 2>&1 "$@"
        }

        pash_redir_all_output_always_execute()
        {
            >>"$PASH_REDIR" 2>&1 "$@"
        }
    fi
fi

export -f pash_redir_output
export -f pash_redir_all_output
export -f pash_redir_all_output_always_execute

source "$PASH_TOP/compiler/orchestrator_runtime/pash_orch_lib.sh"


if [ "$pash_daemon_communicates_through_unix_pipes_flag" -eq 1 ]; then
    pash_communicate_daemon()
    {
        local message=$1
        pash_redir_output echo "Sending msg to daemon: $message"
        echo "$message" > "$RUNTIME_IN_FIFO"
        daemon_response=$(cat "$RUNTIME_OUT_FIFO")
        pash_redir_output echo "Got response from daemon: $daemon_response"
        echo "$daemon_response"
    }

    pash_communicate_daemon_just_send()
    {
        local message=$1
        pash_redir_output echo "Sending msg to daemon: $message"
        echo "$message" > "$RUNTIME_IN_FIFO"
    }

    pash_wait_until_daemon_listening()
    {
        :
    }
else
    pash_communicate_daemon()
    {
        local message=$1
        pash_communicate_unix_socket "compilation-server" "${DAEMON_SOCKET}" "${message}"
    }

    pash_communicate_daemon_just_send()
    {
        pash_communicate_daemon "$1"
    }

    pash_wait_until_daemon_listening()
    {
        pash_wait_until_unix_socket_listening "compilation-server" "${DAEMON_SOCKET}"
    }
fi
export -f pash_communicate_daemon
export -f pash_communicate_daemon_just_send
export -f pash_wait_until_daemon_listening


if [ "$distributed_exec" -eq 1 ]; then
    pash_communicate_worker_manager()
    {
        local message=$1
        pash_redir_output echo "Sending msg to worker manager: $message"
        manager_response=$(echo "$message" | nc -U "$DSPASH_SOCKET")
        pash_redir_output echo "Got response from worker manager: $manager_response"
        echo "$manager_response"
    }
fi

if [ "${pash_speculative_flag}" -eq 1 ]; then
    source "$RUNTIME_DIR/speculative/pash_spec_init_setup.sh"
else
    ## Normal PaSh mode
    ## Exports $daemon_pid
    start_server()
    {
        python3 -S "$PASH_TOP/compiler/pash_compilation_server.py" "$@" &
        export daemon_pid=$!
        ## Wait until daemon has established connection
        pash_wait_until_daemon_listening
    }

    cleanup_server()
    {
        local daemon_pid=$1
        ## Only wait for daemon if it lives (it might be dead, rip)
        if ps -p "$daemon_pid" > /dev/null
        then
            ## Send and receive from daemon
            msg="Done"
            daemon_response=$(pash_communicate_daemon "$msg")
            if [ "$distributed_exec" -eq 1 ]; then
                # kill $worker_manager_pid
                manager_response=$(pash_communicate_worker_manager "$msg")
            fi
            wait 2> /dev/null 1>&2 
        fi
    }
fi

