export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_wHLNLXG/ 
mkdir -p /tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/ 
mkdir -p /tmp/pash_wHLNLXG/056660446ced4e41ab319536b82b8ade/ 
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
{ rm -f "/tmp/pash_wHLNLXG/056660446ced4e41ab319536b82b8ade/#fifo34" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/056660446ced4e41ab319536b82b8ade/#fifo50" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/056660446ced4e41ab319536b82b8ade/#fifo66" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo136" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo137" ; } 
 { rm -f "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo198" ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_wHLNLXG/056660446ced4e41ab319536b82b8ade/#fifo34" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/056660446ced4e41ab319536b82b8ade/#fifo50" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/056660446ced4e41ab319536b82b8ade/#fifo66" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo136" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo137" ; } 
 { mkfifo "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo198" ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_wrap bash -c ' cut -c 88-92 ' <"/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo137" >"/tmp/pash_wHLNLXG/056660446ced4e41ab319536b82b8ade/#fifo34" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_wrap bash -c ' grep -v 999 ' <"/tmp/pash_wHLNLXG/056660446ced4e41ab319536b82b8ade/#fifo34" >"/tmp/pash_wHLNLXG/056660446ced4e41ab319536b82b8ade/#fifo50" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_unwrap <"/tmp/pash_wHLNLXG/056660446ced4e41ab319536b82b8ade/#fifo50" >"/tmp/pash_wHLNLXG/056660446ced4e41ab319536b82b8ade/#fifo66" & }
pids_to_kill="${!} ${pids_to_kill}"
{ sort -r -n <"/tmp/pash_wHLNLXG/056660446ced4e41ab319536b82b8ade/#fifo66" >"/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo198" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo136" -o "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo137" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*d72d9108-355b-4afd-9dd3-2b1bbadd55f9*1*0*/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo136 send*1d82d0ef-3f2e-43cf-881d-207a829976de*0*1*/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo198 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
