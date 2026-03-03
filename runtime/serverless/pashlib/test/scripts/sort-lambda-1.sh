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
{ python3 aws/s3-chunk-reader-approx-correction.py "oneliners/inputs/100M.txt" "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo18" '[{"start": 0, "end": 3276802, "block_id": 0, "shard_id": 0}, {"start": 6553606, "end": 9830408, "block_id": 2, "shard_id": 1}, {"start": 13107212, "end": 16384014, "block_id": 4, "shard_id": 2}, {"start": 19660818, "end": 22937620, "block_id": 6, "shard_id": 3}, {"start": 26214424, "end": 29491226, "block_id": 8, "shard_id": 4}, {"start": 32768030, "end": 36044832, "block_id": 10, "shard_id": 5}, {"start": 39321636, "end": 42598438, "block_id": 12, "shard_id": 6}, {"start": 45875242, "end": 49152044, "block_id": 14, "shard_id": 7}, {"start": 52428848, "end": 55705650, "block_id": 16, "shard_id": 8}, {"start": 58982454, "end": 62259256, "block_id": 18, "shard_id": 9}, {"start": 65536060, "end": 68812862, "block_id": 20, "shard_id": 10}, {"start": 72089666, "end": 75366468, "block_id": 22, "shard_id": 11}, {"start": 78643272, "end": 81920074, "block_id": 24, "shard_id": 12}, {"start": 85196878, "end": 88473680, "block_id": 26, "shard_id": 13}, {"start": 91750484, "end": 95027286, "block_id": 28, "shard_id": 14}, {"start": 98304090, "end": 101580892, "block_id": 30, "shard_id": 15}]' shard=0 num_shards=2 job_uid=74e4fe8f-f4ab-4c17-99a8-6280c4889666 debug=True window_size=None chunks_per_lambda=16 write_headers=false>"/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo18" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/pashlib-ft send*12329258-0e5c-4b3f-a37e-613d187b9961*0*1*/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo9 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )