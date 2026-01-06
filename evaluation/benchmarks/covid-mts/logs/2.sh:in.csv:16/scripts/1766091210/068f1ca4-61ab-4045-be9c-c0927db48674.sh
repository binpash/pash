export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_cahH684/ 
mkdir -p /tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/ 
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
{ rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo433" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo434" ; } 
 { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo495" ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo433" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo434" ; } 
 { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo495" ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ sort -k 1 -n <"/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo434" >"/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo495" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo433" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo434" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*9af4dae8-3f34-481e-ab00-f5dcc71c75a4*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo433 send*aacd2f68-7a0b-40fc-abee-a226153de89f*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo495 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
