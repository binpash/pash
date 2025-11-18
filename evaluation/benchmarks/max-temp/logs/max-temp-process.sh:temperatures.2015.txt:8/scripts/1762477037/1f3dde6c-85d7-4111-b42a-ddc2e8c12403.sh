export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_1hm9dlJ/ 
mkdir -p /tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/ 
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
{ rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo116" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo117" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo120" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo121" ; } 
 { rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo134" ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo116" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo117" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo120" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo121" ; } 
 { mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo134" ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ sort -r -n -m "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo117" "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo121" >"/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo134" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo116" -o "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo117" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo120" -o "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo121" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*9abeab00-9ae4-4121-b6a7-bf1e30d36716*1*0*/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo116 recv*0711c13f-3aae-4d16-9249-774746a29994*1*0*/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo120 send*7096b33b-0d6f-44f0-8f53-254ba34fd86d*0*1*/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo134 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
