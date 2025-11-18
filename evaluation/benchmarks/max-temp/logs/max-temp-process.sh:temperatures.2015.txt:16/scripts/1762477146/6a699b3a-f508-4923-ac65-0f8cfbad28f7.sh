export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_wHLNLXG/ 
mkdir -p /tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/ 
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

rm_pash_fifos() {
{ rm -f "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo228" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo229" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo232" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo233" ; } 
 { rm -f "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo262" ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo228" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo229" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo232" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo233" ; } 
 { mkfifo "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo262" ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ sort -r -n -m "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo229" "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo233" >"/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo262" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo228" -o "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo229" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo232" -o "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo233" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*d03fcb03-963b-4f9b-afe3-f8575f8986e1*1*0*/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo228 recv*0e8c1b69-cf2d-488e-89a9-85a9f6fc90be*1*0*/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo232 send*07acd4cc-4979-4a28-9ddf-cb6a898934be*0*1*/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo262 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
