#!/usr/bin/env bash

export PASH_TOP=${PASH_TOP:-${BASH_SOURCE%/*}}
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/lib/"

## Register the signal handlers
trap kill_all SIGTERM SIGINT

kill_all() {
    kill -s SIGKILL 0
    kill -s SIGKILL "$daemon_pid"
}

## Save the umask to first create some files and then revert it
old_umask=$(umask)
umask u=rwx,g=rx,o=rx

## Get bash version for pash
export PASH_BASH_VERSION="${BASH_VERSINFO[@]:0:3}"

## Create temporary directory and communication FIFOs
export PASH_TMP_PREFIX="$(mktemp -d /tmp/pash_XXXXXXX)/"
export PASH_TIMESTAMP="$(date +"%y-%m-%d-%T")"
export RUNTIME_IN_FIFO="${PASH_TMP_PREFIX}/runtime_in_fifo"
export RUNTIME_OUT_FIFO="${PASH_TMP_PREFIX}/runtime_out_fifo"
rm -f "$RUNTIME_IN_FIFO" "$RUNTIME_OUT_FIFO"
mkfifo "$RUNTIME_IN_FIFO" "$RUNTIME_OUT_FIFO"
export DAEMON_SOCKET="${PASH_TMP_PREFIX}/daemon_socket"
export DSPASH_SOCKET="${PASH_TMP_PREFIX}/dspash_socket"

###############################################################################
# Argument Parsing
###############################################################################

## Initialize all argument variables with default values
pash_init_arg_defaults() {
    # Script execution args
    input_script=""
    shell_name="pash"
    declare -g -a script_args=()
    command_mode=""
    command_text=""

    # Shell flags
    allexport_flag="+a"
    verbose_flag=""
    xtrace_flag=""

    # Shared arguments (preprocessor + server)
    arg_debug=""
    arg_log_file=""
    arg_speculative=""
    arg_config_path=""
    arg_local_annotations_dir=""

    # Preprocessor-specific arguments
    arg_bash=""

    # Server-specific arguments
    arg_width=""
    arg_no_optimize=""
    arg_dry_run_compiler=""
    arg_assert_compiler_success=""
    arg_assert_all_regions_parallelizable=""
    arg_avoid_pash_runtime_completion=""
    arg_output_optimized=""
    arg_graphviz=""
    arg_graphviz_dir=""
    arg_no_parallel_pipelines=""
    arg_parallel_pipelines_limit=""
    arg_r_split_batch_size=""
    arg_version=""
    arg_no_eager=""
    arg_profile_driven=""
    arg_termination=""
    arg_daemon_communicates_through_unix_pipes=""
    arg_distributed_exec=""

    # Obsolete args (still supported for backward compatibility)
    arg_no_daemon=""
    arg_parallel_pipelines=""
    arg_r_split=""
    arg_dgsh_tee=""
    arg_speculation=""
}

## Parse all command-line arguments into variables
pash_parse_args() {
    local i=1
    local arg next_i next_arg

    while [ $i -le $# ]; do
        arg="${!i}"
        next_i=$((i+1))
        next_arg="${!next_i}"

        case "$arg" in
            # Shell flags
            -c|--command)
                command_mode="-c"
                command_text="$next_arg"
                i=$next_i
                ;;
            -a)
                allexport_flag="-a"
                ;;
            +a)
                allexport_flag="+a"
                ;;
            -v)
                verbose_flag="-v"
                ;;
            -x)
                xtrace_flag="-x"
                ;;

            # Shared arguments (preprocessor + server)
            -d|--debug)
                arg_debug="$next_arg"
                i=$next_i
                ;;
            --log_file)
                arg_log_file="$next_arg"
                i=$next_i
                ;;
            --speculative)
                arg_speculative="--speculative"
                ;;
            --config_path)
                arg_config_path="$next_arg"
                i=$next_i
                ;;
            --local-annotations-dir)
                arg_local_annotations_dir="$next_arg"
                i=$next_i
                ;;

            # Preprocessor-specific
            --bash)
                arg_bash="--bash"
                ;;

            # Server-specific arguments
            -w|--width)
                arg_width="$next_arg"
                i=$next_i
                ;;
            --no_optimize)
                arg_no_optimize="--no_optimize"
                ;;
            --dry_run_compiler)
                arg_dry_run_compiler="--dry_run_compiler"
                ;;
            --assert_compiler_success)
                arg_assert_compiler_success="--assert_compiler_success"
                ;;
            --assert_all_regions_parallelizable)
                arg_assert_all_regions_parallelizable="--assert_all_regions_parallelizable"
                ;;
            --avoid_pash_runtime_completion)
                arg_avoid_pash_runtime_completion="--avoid_pash_runtime_completion"
                ;;
            -p|--output_optimized)
                arg_output_optimized="--output_optimized"
                ;;
            --graphviz)
                arg_graphviz="$next_arg"
                i=$next_i
                ;;
            --graphviz_dir)
                arg_graphviz_dir="$next_arg"
                i=$next_i
                ;;
            --no_parallel_pipelines)
                arg_no_parallel_pipelines="--no_parallel_pipelines"
                ;;
            --parallel_pipelines_limit)
                arg_parallel_pipelines_limit="$next_arg"
                i=$next_i
                ;;
            --r_split_batch_size)
                arg_r_split_batch_size="$next_arg"
                i=$next_i
                ;;
            --version)
                arg_version="--version"
                ;;
            --no_eager)
                arg_no_eager="--no_eager"
                ;;
            --profile_driven)
                arg_profile_driven="--profile_driven"
                ;;
            --termination)
                arg_termination="$next_arg"
                i=$next_i
                ;;
            --daemon_communicates_through_unix_pipes)
                arg_daemon_communicates_through_unix_pipes="--daemon_communicates_through_unix_pipes"
                ;;
            --distributed_exec)
                arg_distributed_exec="--distributed_exec"
                ;;

            # Obsolete arguments (still accept for backward compatibility)
            --no_daemon)
                arg_no_daemon="--no_daemon"
                ;;
            --parallel_pipelines)
                arg_parallel_pipelines="--parallel_pipelines"
                ;;
            --r_split)
                arg_r_split="--r_split"
                ;;
            --dgsh_tee)
                arg_dgsh_tee="--dgsh_tee"
                ;;
            --speculation)
                arg_speculation="$next_arg"
                i=$next_i
                ;;
            -t|--output_time)
                # Obsolete, time is always logged now
                ;;

            # Positional arguments (input script and script args)
            *)
                if [ -z "$input_script" ] && [ -z "$command_mode" ]; then
                    input_script="$arg"
                    shell_name="$arg"
                else
                    script_args+=("$arg")
                fi
                ;;
        esac
        i=$((i+1))
    done

    # In -c mode, first positional arg becomes shell_name ($0)
    if [ -n "$command_mode" ] && [ ${#script_args[@]} -gt 0 ]; then
        shell_name="${script_args[0]}"
        script_args=("${script_args[@]:1}")
    fi
}

## Handle -c command mode: write command text to a temp file
pash_handle_command_mode() {
    if [ -n "$command_text" ]; then
        local command_file
        command_file=$(mktemp "${PASH_TMP_PREFIX}/command_XXXXXX.sh")
        printf '%s' "$command_text" > "$command_file"
        input_script="$command_file"
    fi
}

## Build the preprocessor arguments array
pash_build_preprocessor_args() {
    declare -g -a preprocessor_args=()
    preprocessor_args+=("--output" "$preprocessed_output")
    [ -n "$arg_debug" ] && preprocessor_args+=("-d" "$arg_debug")
    [ -n "$arg_log_file" ] && preprocessor_args+=("--log_file" "$arg_log_file")
    [ -n "$arg_bash" ] && preprocessor_args+=("$arg_bash")
    [ -n "$arg_speculative" ] && preprocessor_args+=("$arg_speculative")
    preprocessor_args+=("$input_script")
}

## Build the compilation server arguments array
pash_build_server_args() {
    declare -g -a server_args=()

    # Shared arguments
    [ -n "$arg_debug" ] && server_args+=("-d" "$arg_debug")
    [ -n "$arg_log_file" ] && server_args+=("--log_file" "$arg_log_file")
    [ -n "$arg_speculative" ] && server_args+=("$arg_speculative")
    [ -n "$arg_config_path" ] && server_args+=("--config_path" "$arg_config_path")

    # Server-specific arguments
    [ -n "$arg_width" ] && server_args+=("-w" "$arg_width")
    [ -n "$arg_no_optimize" ] && server_args+=("$arg_no_optimize")
    [ -n "$arg_dry_run_compiler" ] && server_args+=("$arg_dry_run_compiler")
    [ -n "$arg_assert_compiler_success" ] && server_args+=("$arg_assert_compiler_success")
    [ -n "$arg_assert_all_regions_parallelizable" ] && server_args+=("$arg_assert_all_regions_parallelizable")
    [ -n "$arg_avoid_pash_runtime_completion" ] && server_args+=("$arg_avoid_pash_runtime_completion")
    [ -n "$arg_output_optimized" ] && server_args+=("$arg_output_optimized")
    [ -n "$arg_graphviz" ] && server_args+=("--graphviz" "$arg_graphviz")
    [ -n "$arg_graphviz_dir" ] && server_args+=("--graphviz_dir" "$arg_graphviz_dir")
    [ -n "$arg_no_parallel_pipelines" ] && server_args+=("$arg_no_parallel_pipelines")
    [ -n "$arg_parallel_pipelines_limit" ] && server_args+=("--parallel_pipelines_limit" "$arg_parallel_pipelines_limit")
    [ -n "$arg_r_split_batch_size" ] && server_args+=("--r_split_batch_size" "$arg_r_split_batch_size")
    [ -n "$arg_version" ] && server_args+=("$arg_version")
    [ -n "$arg_no_eager" ] && server_args+=("$arg_no_eager")
    [ -n "$arg_profile_driven" ] && server_args+=("$arg_profile_driven")
    [ -n "$arg_termination" ] && server_args+=("--termination" "$arg_termination")
    [ -n "$arg_daemon_communicates_through_unix_pipes" ] && server_args+=("$arg_daemon_communicates_through_unix_pipes")
    [ -n "$arg_distributed_exec" ] && server_args+=("$arg_distributed_exec")

    # Obsolete arguments
    [ -n "$arg_no_daemon" ] && server_args+=("$arg_no_daemon")
    [ -n "$arg_parallel_pipelines" ] && server_args+=("$arg_parallel_pipelines")
    [ -n "$arg_r_split" ] && server_args+=("$arg_r_split")
    [ -n "$arg_dgsh_tee" ] && server_args+=("$arg_dgsh_tee")
    [ -n "$arg_speculation" ] && server_args+=("--speculation" "$arg_speculation")
}

###############################################################################
# Environment Setup (previously in pash_init_setup.sh)
###############################################################################

pash_setup_environment() {
    # Source local pash config if it exists
    ## TODO: Is this used ever?
    [ -f ~/.pash_init ] && source ~/.pash_init

    # Directory exports
    export RUNTIME_DIR="$PASH_TOP/compiler/orchestrator_runtime"
    export WRAPPER_LIB_DIR="$RUNTIME_DIR/../wrapper_library/"
    export RUNTIME_LIBRARY_DIR="$RUNTIME_DIR/../../runtime/"

    # Debug/logging setup
    export PASH_DEBUG_LEVEL="${arg_debug:-0}"
    if [ -n "$arg_log_file" ]; then
        export PASH_REDIR="$arg_log_file"
    else
        export PASH_REDIR="&2"
    fi

    # Local annotations setup
    if [ -n "$arg_local_annotations_dir" ]; then
        export ANNOTATIONS_PATH=$(realpath "$arg_local_annotations_dir")
        export PYTHONPATH="$ANNOTATIONS_PATH:$PYTHONPATH"
    fi

    # Export flags used by runtime scripts
    export pash_output_time_flag=1
    export pash_execute_flag=1
    export pash_dry_run_compiler_flag=$( [ -n "$arg_dry_run_compiler" ] && echo 1 || echo 0 )
    export pash_assert_compiler_success_flag=$( [ -n "$arg_assert_compiler_success" ] && echo 1 || echo 0 )
    export pash_assert_all_regions_parallelizable_flag=$( [ -n "$arg_assert_all_regions_parallelizable" ] && echo 1 || echo 0 )
    export pash_avoid_pash_runtime_completion_flag=$( [ -n "$arg_avoid_pash_runtime_completion" ] && echo 1 || echo 0 )
    export pash_profile_driven_flag=$( [ -n "$arg_profile_driven" ] && echo 1 || echo 1 )  # Default 1
    export pash_no_parallel_pipelines=$( [ -n "$arg_no_parallel_pipelines" ] && echo 1 || echo 0 )
    export pash_daemon_communicates_through_unix_pipes_flag=$( [ -n "$arg_daemon_communicates_through_unix_pipes" ] && echo 1 || echo 0 )
    export pash_speculative_flag=$( [ -n "$arg_speculative" ] && echo 1 || echo 0 )
    export show_version=$( [ -n "$arg_version" ] && echo 1 || echo 0 )
    export distributed_exec=$( [ -n "$arg_distributed_exec" ] && echo 1 || echo 0 )
    [ -n "$arg_bash" ] && export pash_shell="bash"
}

###############################################################################
# Logging Functions (previously in pash_init_setup.sh)
###############################################################################

pash_setup_logging() {
    if [ "$PASH_DEBUG_LEVEL" -eq 0 ]; then
        pash_redir_output() {
            :
        }

        pash_redir_all_output() {
            :
        }

        pash_redir_all_output_always_execute() {
            > /dev/null 2>&1 "$@"
        }
    else
        if [ "$PASH_REDIR" == '&2' ]; then
            pash_redir_output() {
                >&2 "$@"
            }

            pash_redir_all_output() {
                >&2 "$@"
            }

            pash_redir_all_output_always_execute() {
                >&2 "$@"
            }
        else
            pash_redir_output() {
                >>"$PASH_REDIR" "$@"
            }

            pash_redir_all_output() {
                >>"$PASH_REDIR" 2>&1 "$@"
            }

            pash_redir_all_output_always_execute() {
                >>"$PASH_REDIR" 2>&1 "$@"
            }
        fi
    fi

    export -f pash_redir_output
    export -f pash_redir_all_output
    export -f pash_redir_all_output_always_execute
}

###############################################################################
# Communication Functions (previously in pash_orch_lib.sh and pash_init_setup.sh)
###############################################################################

pash_setup_communication() {
    ## Socket communication utilities (from pash_orch_lib.sh)
    pash_wait_until_unix_socket_listening() {
        local server_name=$1
        local socket=$2
        local i=0
        local maximum_retries=1000
        until echo "Daemon Start" 2> /dev/null | nc -U "$socket" >/dev/null 2>&1; do
            sleep 0.01
            i=$((i+1))
            if [ $i -eq $maximum_retries ]; then
                echo "Error: Maximum retries: $maximum_retries exceeded when waiting for server: ${server_name} to bind to socket: ${socket}!" 1>&2
                echo "Exiting..." 1>&2
                exit 1
            fi
        done
    }

    pash_communicate_unix_socket() {
        local server_name=$1
        local socket=$2
        local message=$3
        pash_redir_output echo "Sending msg to ${server_name}: $message"
        daemon_response=$(echo "$message" | nc -U "${socket}")
        pash_redir_output echo "Got response from ${server_name}: $daemon_response"
        echo "$daemon_response"
    }

    export -f pash_wait_until_unix_socket_listening
    export -f pash_communicate_unix_socket

    ## Daemon communication (conditional on unix pipes flag)
    if [ "$pash_daemon_communicates_through_unix_pipes_flag" -eq 1 ]; then
        pash_communicate_daemon() {
            local message=$1
            pash_redir_output echo "Sending msg to daemon: $message"
            echo "$message" > "$RUNTIME_IN_FIFO"
            daemon_response=$(cat "$RUNTIME_OUT_FIFO")
            pash_redir_output echo "Got response from daemon: $daemon_response"
            echo "$daemon_response"
        }

        pash_communicate_daemon_just_send() {
            local message=$1
            pash_redir_output echo "Sending msg to daemon: $message"
            echo "$message" > "$RUNTIME_IN_FIFO"
        }

        pash_wait_until_daemon_listening() {
            :
        }
    else
        pash_communicate_daemon() {
            local message=$1
            pash_communicate_unix_socket "compilation-server" "${DAEMON_SOCKET}" "${message}"
        }

        pash_communicate_daemon_just_send() {
            pash_communicate_daemon "$1"
        }

        pash_wait_until_daemon_listening() {
            pash_wait_until_unix_socket_listening "compilation-server" "${DAEMON_SOCKET}"
        }
    fi

    export -f pash_communicate_daemon
    export -f pash_communicate_daemon_just_send
    export -f pash_wait_until_daemon_listening

    ## Worker manager communication (for distributed execution)
    if [ "$distributed_exec" -eq 1 ]; then
        pash_communicate_worker_manager() {
            local message=$1
            pash_redir_output echo "Sending msg to worker manager: $message"
            manager_response=$(echo "$message" | nc -U "$DSPASH_SOCKET")
            pash_redir_output echo "Got response from worker manager: $manager_response"
            echo "$manager_response"
        }
        export -f pash_communicate_worker_manager
    fi
}

###############################################################################
# Server Functions
###############################################################################

pash_setup_server_functions() {
    if [ "${pash_speculative_flag}" -eq 1 ]; then
        ## Speculative mode - source the speculative setup
        source "$PASH_TOP/compiler/orchestrator_runtime/speculative/pash_spec_init_setup.sh"
    else
        ## Normal PaSh mode
        start_server() {
            "$PASH_TOP/python_pkgs/bin/python" "$PASH_TOP/compiler/pash_compilation_server.py" "$@" &
            export daemon_pid=$!
            pash_wait_until_daemon_listening
        }

        cleanup_server() {
            local daemon_pid=$1
            if ps -p "$daemon_pid" > /dev/null 2>&1; then
                msg="Done"
                daemon_response=$(pash_communicate_daemon "$msg")
                if [ "$distributed_exec" -eq 1 ]; then
                    manager_response=$(pash_communicate_worker_manager "$msg")
                fi
                wait 2> /dev/null 1>&2
            fi
        }

        export -f start_server
        export -f cleanup_server
    fi
}

###############################################################################
# Main Execution
###############################################################################

## 1. Parse arguments
pash_init_arg_defaults
pash_parse_args "$@"

## 2. Setup environment and functions based on parsed arguments
pash_setup_environment
pash_setup_logging
pash_setup_communication
pash_setup_server_functions

## 3. Handle -c command mode
pash_handle_command_mode

## 4. Create temporary file for preprocessed output
preprocessed_output=$(mktemp "${PASH_TMP_PREFIX}/preprocessed_XXXXXX.sh")

## 5. Build argument arrays for preprocessor and server
pash_build_preprocessor_args
pash_build_server_args

## 6. Start the compilation server
if [ "$show_version" -eq 0 ]; then
    start_server "${server_args[@]}"
fi

## 7. Restore umask before executing user scripts
umask "$old_umask"

## 8. Run the preprocessor
PASH_FROM_SH="PaSh preprocessor" "$PASH_TOP/python_pkgs/bin/python" \
    "$PASH_TOP/compiler/pash_preprocessor.py" "${preprocessor_args[@]}"
pash_exit_code=$?

## 9. If preprocessing succeeded, execute the preprocessed script
if [ "$pash_exit_code" -eq 0 ]; then
    bash_flags="$allexport_flag $verbose_flag $xtrace_flag"
    # shellcheck disable=SC2086
    bash $bash_flags -c "source $preprocessed_output" "$shell_name" "${script_args[@]}"
    pash_exit_code=$?
fi

## 10. Cleanup
rm -f "$preprocessed_output"
if [ "$show_version" -eq 0 ]; then
    cleanup_server "${daemon_pid}"
fi

if [ "$PASH_DEBUG_LEVEL" -eq 0 ]; then
    rm -rf "${PASH_TMP_PREFIX}"
fi

(exit "$pash_exit_code")
