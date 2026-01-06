export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_ybt07Kl/ 
mkdir -p /tmp/pash_ybt07Kl/ac90e51ec3f54dd5b02e041d07d42628/ 
mkdir -p /tmp/pash_ybt07Kl/e9be7e5d448d4dad9e4d8eee8653328a/ 
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
{ rm -f "/tmp/pash_ybt07Kl/e9be7e5d448d4dad9e4d8eee8653328a/#fifo6" ; } 
 { { rm -f "/tmp/pash_ybt07Kl/e9be7e5d448d4dad9e4d8eee8653328a/#fifo8" ; } 
 { { rm -f "/tmp/pash_ybt07Kl/ac90e51ec3f54dd5b02e041d07d42628/#fifo33" ; } 
 { { rm -f "/tmp/pash_ybt07Kl/ac90e51ec3f54dd5b02e041d07d42628/#fifo34" ; } 
 { { rm -f "/tmp/pash_ybt07Kl/ac90e51ec3f54dd5b02e041d07d42628/#fifo37" ; } 
 { { rm -f "/tmp/pash_ybt07Kl/ac90e51ec3f54dd5b02e041d07d42628/#fifo38" ; } 
 { { rm -f "/tmp/pash_ybt07Kl/ac90e51ec3f54dd5b02e041d07d42628/#fifo39" ; } 
 { rm -f "/tmp/pash_ybt07Kl/ac90e51ec3f54dd5b02e041d07d42628/#fifo43" ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_ybt07Kl/e9be7e5d448d4dad9e4d8eee8653328a/#fifo6" ; } 
 { { mkfifo "/tmp/pash_ybt07Kl/e9be7e5d448d4dad9e4d8eee8653328a/#fifo8" ; } 
 { { mkfifo "/tmp/pash_ybt07Kl/ac90e51ec3f54dd5b02e041d07d42628/#fifo33" ; } 
 { { mkfifo "/tmp/pash_ybt07Kl/ac90e51ec3f54dd5b02e041d07d42628/#fifo34" ; } 
 { { mkfifo "/tmp/pash_ybt07Kl/ac90e51ec3f54dd5b02e041d07d42628/#fifo37" ; } 
 { { mkfifo "/tmp/pash_ybt07Kl/ac90e51ec3f54dd5b02e041d07d42628/#fifo38" ; } 
 { { mkfifo "/tmp/pash_ybt07Kl/ac90e51ec3f54dd5b02e041d07d42628/#fifo39" ; } 
 { mkfifo "/tmp/pash_ybt07Kl/ac90e51ec3f54dd5b02e041d07d42628/#fifo43" ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_merge "/tmp/pash_ybt07Kl/ac90e51ec3f54dd5b02e041d07d42628/#fifo34" "/tmp/pash_ybt07Kl/ac90e51ec3f54dd5b02e041d07d42628/#fifo38" >"/tmp/pash_ybt07Kl/e9be7e5d448d4dad9e4d8eee8653328a/#fifo6" & }
pids_to_kill="${!} ${pids_to_kill}"
{ sort -u <"/tmp/pash_ybt07Kl/e9be7e5d448d4dad9e4d8eee8653328a/#fifo6" >"/tmp/pash_ybt07Kl/e9be7e5d448d4dad9e4d8eee8653328a/#fifo8" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_split "/tmp/pash_ybt07Kl/e9be7e5d448d4dad9e4d8eee8653328a/#fifo8" 1000000 "/tmp/pash_ybt07Kl/ac90e51ec3f54dd5b02e041d07d42628/#fifo39" "/tmp/pash_ybt07Kl/ac90e51ec3f54dd5b02e041d07d42628/#fifo43" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_ybt07Kl/ac90e51ec3f54dd5b02e041d07d42628/#fifo33" -o "/tmp/pash_ybt07Kl/ac90e51ec3f54dd5b02e041d07d42628/#fifo34" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_ybt07Kl/ac90e51ec3f54dd5b02e041d07d42628/#fifo37" -o "/tmp/pash_ybt07Kl/ac90e51ec3f54dd5b02e041d07d42628/#fifo38" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*84437d8a-6122-4f61-8d52-debc38b458e5*1*0*/tmp/pash_ybt07Kl/ac90e51ec3f54dd5b02e041d07d42628/#fifo33 recv*2c321ff9-0154-42b2-85bc-28b1172429a8*1*0*/tmp/pash_ybt07Kl/ac90e51ec3f54dd5b02e041d07d42628/#fifo37 send*c4a5c0a8-1a38-4b88-9b30-d83f814ab94a*0*1*/tmp/pash_ybt07Kl/ac90e51ec3f54dd5b02e041d07d42628/#fifo39 send*683511aa-2ff1-4bfc-8a87-7dafa2bb1e8a*0*1*/tmp/pash_ybt07Kl/ac90e51ec3f54dd5b02e041d07d42628/#fifo43 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
