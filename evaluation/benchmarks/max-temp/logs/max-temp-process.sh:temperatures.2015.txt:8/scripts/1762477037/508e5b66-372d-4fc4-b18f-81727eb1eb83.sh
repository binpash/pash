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
{ rm -f "/tmp/pash_1hm9dlJ/9f268ab6ed4241d59b51afd41f5c68b9/#fifo22" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/9f268ab6ed4241d59b51afd41f5c68b9/#fifo30" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/9f268ab6ed4241d59b51afd41f5c68b9/#fifo38" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo72" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo73" ; } 
 { rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo102" ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_1hm9dlJ/9f268ab6ed4241d59b51afd41f5c68b9/#fifo22" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/9f268ab6ed4241d59b51afd41f5c68b9/#fifo30" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/9f268ab6ed4241d59b51afd41f5c68b9/#fifo38" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo72" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo73" ; } 
 { mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo102" ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_wrap bash -c ' cut -c 88-92 ' <"/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo73" >"/tmp/pash_1hm9dlJ/9f268ab6ed4241d59b51afd41f5c68b9/#fifo22" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_wrap bash -c ' grep -v 999 ' <"/tmp/pash_1hm9dlJ/9f268ab6ed4241d59b51afd41f5c68b9/#fifo22" >"/tmp/pash_1hm9dlJ/9f268ab6ed4241d59b51afd41f5c68b9/#fifo30" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_unwrap <"/tmp/pash_1hm9dlJ/9f268ab6ed4241d59b51afd41f5c68b9/#fifo30" >"/tmp/pash_1hm9dlJ/9f268ab6ed4241d59b51afd41f5c68b9/#fifo38" & }
pids_to_kill="${!} ${pids_to_kill}"
{ sort -r -n <"/tmp/pash_1hm9dlJ/9f268ab6ed4241d59b51afd41f5c68b9/#fifo38" >"/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo102" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo72" -o "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo73" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*6d7b68bd-6710-4413-9787-30c19de963b1*1*0*/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo72 send*b597cdd2-ecf6-4bf3-bc86-6c2fa5651fa3*0*1*/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo102 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
