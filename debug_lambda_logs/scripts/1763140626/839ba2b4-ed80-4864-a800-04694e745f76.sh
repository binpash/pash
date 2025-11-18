export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_p4qiBDD/ 
mkdir -p /tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/ 
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
{ rm -f "/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo20" ; } 
 { { rm -f "/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo21" ; } 
 { { rm -f "/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo22" ; } 
 { { rm -f "/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo27" ; } 
 { rm -f "/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo28" ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo20" ; } 
 { { mkfifo "/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo21" ; } 
 { { mkfifo "/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo22" ; } 
 { { mkfifo "/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo27" ; } 
 { mkfifo "/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo28" ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ sort -m "/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo21" "/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo28" >"/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo22" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo20" -o "/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo21" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo27" -o "/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo28" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*f0082947-1beb-4dda-b5f4-3e82e8a2ef0b*1*0*/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo20 recv*61d430ca-72a0-4745-bf97-94441b9af0fc*1*0*/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo27 & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3.9 aws/s3-put-object.py oneliners/outputs/sort.sh:1M.txt:2:leashstdout.txt "/tmp/pash_p4qiBDD/5aeeab8868424403b26ff7b77f815013/#fifo22" $1 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
