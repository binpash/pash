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
{ rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo108" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo109" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo112" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo113" ; } 
 { rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo130" ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo108" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo109" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo112" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo113" ; } 
 { mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo130" ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ sort -r -n -m "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo109" "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo113" >"/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo130" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo108" -o "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo109" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo112" -o "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo113" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*26fe3ace-609e-4efc-aa02-513f57bc73d4*1*0*/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo108 recv*93c44d69-de53-4147-8460-2eb7e12dbb99*1*0*/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo112 send*269681cb-aa28-4ca9-8dab-ceb8a793729d*0*1*/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo130 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
