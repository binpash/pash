export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_vRTnkmJ/ 
mkdir -p /tmp/pash_vRTnkmJ/98566d95cd3246499e2519ef123faa17/ 
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
{ rm -f "/tmp/pash_vRTnkmJ/98566d95cd3246499e2519ef123faa17/#fifo18" ; } 
 { { rm -f "/tmp/pash_vRTnkmJ/98566d95cd3246499e2519ef123faa17/#fifo19" ; } 
 { rm -f "/tmp/pash_vRTnkmJ/98566d95cd3246499e2519ef123faa17/#fifo32" ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_vRTnkmJ/98566d95cd3246499e2519ef123faa17/#fifo18" ; } 
 { { mkfifo "/tmp/pash_vRTnkmJ/98566d95cd3246499e2519ef123faa17/#fifo19" ; } 
 { mkfifo "/tmp/pash_vRTnkmJ/98566d95cd3246499e2519ef123faa17/#fifo32" ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ sort <"/tmp/pash_vRTnkmJ/98566d95cd3246499e2519ef123faa17/#fifo19" >"/tmp/pash_vRTnkmJ/98566d95cd3246499e2519ef123faa17/#fifo32" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vRTnkmJ/98566d95cd3246499e2519ef123faa17/#fifo18" -o "/tmp/pash_vRTnkmJ/98566d95cd3246499e2519ef123faa17/#fifo19" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*c2eb89e5-5bbb-43a8-a47a-32a3d0182676*1*0*/tmp/pash_vRTnkmJ/98566d95cd3246499e2519ef123faa17/#fifo18 send*ba3ec213-fe9b-4849-8246-93447919cfc9*0*1*/tmp/pash_vRTnkmJ/98566d95cd3246499e2519ef123faa17/#fifo32 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
