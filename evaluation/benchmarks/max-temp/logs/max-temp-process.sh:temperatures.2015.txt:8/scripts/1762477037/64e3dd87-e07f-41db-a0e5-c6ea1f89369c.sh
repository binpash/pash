export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_1hm9dlJ/ 
mkdir -p /tmp/pash_1hm9dlJ/9f268ab6ed4241d59b51afd41f5c68b9/ 
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
{ rm -f "/tmp/pash_1hm9dlJ/9f268ab6ed4241d59b51afd41f5c68b9/#fifo26" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/9f268ab6ed4241d59b51afd41f5c68b9/#fifo34" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/9f268ab6ed4241d59b51afd41f5c68b9/#fifo42" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo88" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo89" ; } 
 { rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo118" ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_1hm9dlJ/9f268ab6ed4241d59b51afd41f5c68b9/#fifo26" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/9f268ab6ed4241d59b51afd41f5c68b9/#fifo34" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/9f268ab6ed4241d59b51afd41f5c68b9/#fifo42" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo88" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo89" ; } 
 { mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo118" ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_wrap bash -c ' cut -c 88-92 ' <"/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo89" >"/tmp/pash_1hm9dlJ/9f268ab6ed4241d59b51afd41f5c68b9/#fifo26" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_wrap bash -c ' grep -v 999 ' <"/tmp/pash_1hm9dlJ/9f268ab6ed4241d59b51afd41f5c68b9/#fifo26" >"/tmp/pash_1hm9dlJ/9f268ab6ed4241d59b51afd41f5c68b9/#fifo34" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_unwrap <"/tmp/pash_1hm9dlJ/9f268ab6ed4241d59b51afd41f5c68b9/#fifo34" >"/tmp/pash_1hm9dlJ/9f268ab6ed4241d59b51afd41f5c68b9/#fifo42" & }
pids_to_kill="${!} ${pids_to_kill}"
{ sort -r -n <"/tmp/pash_1hm9dlJ/9f268ab6ed4241d59b51afd41f5c68b9/#fifo42" >"/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo118" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo88" -o "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo89" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*f5a8fd98-2c08-406f-b273-b92f05cd413d*1*0*/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo88 send*0711c13f-3aae-4d16-9249-774746a29994*0*1*/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo118 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
