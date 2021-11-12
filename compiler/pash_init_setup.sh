# source the local pash config
source ~/.pash_init
## File directory
export RUNTIME_DIR=$(dirname "${BASH_SOURCE[0]}")
## TODO: Is there a better way to do this?
export RUNTIME_LIBRARY_DIR="$RUNTIME_DIR/../runtime/"
export PASH_REDIR="&2"
export PASH_DEBUG_LEVEL=0

## Check flags
export pash_output_time_flag=1
export pash_execute_flag=1
export pash_speculation_flag=0 # By default there is no speculation
export pash_dry_run_compiler_flag=0
export pash_assert_compiler_success_flag=0
export pash_checking_speculation=0
export pash_checking_log_file=0
export pash_checking_debug_level=0
export pash_avoid_pash_runtime_completion_flag=0
export pash_parallel_pipelines=0

for item in $@
do
    if [ "$pash_checking_speculation" -eq 1 ]; then
        export pash_checking_speculation=0
        if [ "no_spec" == "$item" ]; then
            export pash_speculation_flag=0
        elif [ "quick_abort" == "$item" ]; then
            ## TODO: Fix how speculation interacts with dry_run, assert_compiler_success
            export pash_speculation_flag=1
        else
            echo "$$: Unknown value for option --speculation" 1>&2
            exit 1
        fi
    fi

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

    if [ "--dry_run_compiler" == "$item" ]; then
        export pash_dry_run_compiler_flag=1
    fi

    if [ "--assert_compiler_success" == "$item" ]; then
        export pash_assert_compiler_success_flag=1
    fi

    if [ "--speculation" == "$item" ]; then
        pash_checking_speculation=1
    fi

    if [ "--log_file" == "$item" ]; then
        pash_checking_log_file=1
    fi

    if [ "--avoid_pash_runtime_completion" == "$item" ]; then
        export pash_avoid_pash_runtime_completion_flag=1
    fi

    if [ "-d" == "$item" ] || [ "--debug" == "$item" ]; then
        pash_checking_debug_level=1
    fi

    ## TODO: Add this flag in pash.py too so that it is printed in help.
    if [ "--pash_parallel_pipelines" == "$item" ]; then
        export pash_parallel_pipelines=1
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
        > /dev/null 2>&1 $@
    }

else
    if [ "$PASH_REDIR" == '&2' ]; then
        pash_redir_output()
        {
            >&2 $@
        }

        pash_redir_all_output()
        {
            >&2 $@
        }

        pash_redir_all_output_always_execute()
        {
            >&2 $@
        }
    else
        pash_redir_output()
        {
            >>"$PASH_REDIR" $@
        }

        pash_redir_all_output()
        {
            >>"$PASH_REDIR" 2>&1 $@
        }

        pash_redir_all_output_always_execute()
        {
            >>"$PASH_REDIR" 2>&1 $@
        }
    fi
fi

export -f pash_redir_output
export -f pash_redir_all_output
export -f pash_redir_all_output_always_execute


pash_send_to_daemon()
{
    local message="$1"
    echo "$message" | nc -U "$DAEMON_SOCKET"
    # echo "$message" > "$RUNTIME_IN_FIFO"
    pash_redir_output echo "Sent msg to daemon: $message"
}

pash_receive_from_daemon()
{
    daemon_response=$(cat "$RUNTIME_OUT_FIFO")
    pash_redir_output echo "Got response from daemon: $daemon_response"
    echo "$daemon_response"
}

## TODO: Make sure to fix that for both pipes and sockets
pash_communicate_daemon()
{
    local message=$1
    pash_redir_output echo "Sending msg to daemon: $message"
    daemon_response=$(echo "$message" | nc -U "$DAEMON_SOCKET")
    pash_redir_output echo "Got response from daemon: $daemon_response"
    echo "$daemon_response"
}

export -f pash_send_to_daemon
export -f pash_receive_from_daemon
export -f pash_communicate_daemon
