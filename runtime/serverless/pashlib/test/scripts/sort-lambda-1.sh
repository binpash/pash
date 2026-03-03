export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export LOCPATH=/var/task/runtime/locale
export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_mExGzpm/ 
mkdir -p /tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/ 
inform_daemon_exit () 
{ 
    msg="Exit:${process_id}";
    daemon_response=$(pash_communicate_daemon_just_send "$msg")
}
pash_communicate_daemon () 
{ 
    local message=$1;
    pash_communicate_unix_socket "compilation-server" "${DAEMON_SOCKET}" "${message}"
}
declare -fx pash_communicate_daemon
pash_communicate_daemon_just_send () 
{ 
    pash_communicate_daemon "$1"
}
declare -fx pash_communicate_daemon_just_send
pash_communicate_unix_socket () 
{ 
    local server_name=$1;
    local socket=$2;
    local message=$3;
    pash_redir_output echo "Sending msg to ${server_name}: $message";
    daemon_response=$(echo "$message" | nc -U "${socket}");
    pash_redir_output echo "Got response from ${server_name}: $daemon_response";
    echo "$daemon_response"
}
declare -fx pash_communicate_unix_socket
pash_redir_all_output () 
{ 
    :
}
declare -fx pash_redir_all_output
pash_redir_all_output_always_execute () 
{ 
    "$@" > /dev/null 2>&1
}
declare -fx pash_redir_all_output_always_execute
pash_redir_output () 
{ 
    :
}
declare -fx pash_redir_output
pash_wait_until_daemon_listening () 
{ 
    pash_wait_until_unix_socket_listening "compilation-server" "${DAEMON_SOCKET}"
}
declare -fx pash_wait_until_daemon_listening
pash_wait_until_unix_socket_listening () 
{ 
    local server_name=$1;
    local socket=$2;
    i=0;
    maximum_retries=1000;
    until echo "Daemon Start" 2> /dev/null | nc -U "$socket" > /dev/null 2>&1; do
        sleep 0.01;
        i=$((i+1));
        if [ $i -eq $maximum_retries ]; then
            echo "Error: Maximum retries: $maximum_retries exceeded when waiting for server: ${server_name} to bind to socket: ${socket}!" 1>&2;
            echo "Exiting..." 1>&2;
            exit 1;
        fi;
    done
}
declare -fx pash_wait_until_unix_socket_listening
run_parallel () 
{ 
    trap inform_daemon_exit SIGTERM SIGINT EXIT;
    export SCRIPT_TO_EXECUTE="$pash_script_to_execute";
    source "$RUNTIME_DIR/pash_restore_state_and_execute.sh"
}

rm_pash_fifos() {
{ rm -f "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo9" ; } 
 { { rm -f "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo18" ; } 
 { rm -f "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo19" ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo9" ; } 
 { { mkfifo "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo18" ; } 
 { mkfifo "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo19" ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ sort <"/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo19" >"/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo9" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo18" -o "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo19" -I -m 1G -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3 aws/s3-chunk-reader-approx-correction.py "oneliners/inputs/1M.txt" "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo18" '[{"start": 0, "end": 524287, "block_id": 0, "shard_id": 0}]' shard=0 num_shards=2 job_uid=efba3ebf-ca02-4246-a4a0-504288102afc debug=True window_size=None chunks_per_lambda=1 write_headers=false >"/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo18" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/pashlib-ft send*12329258-0e5c-4b3f-a37e-613d187b9961*0*1*/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo9 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )