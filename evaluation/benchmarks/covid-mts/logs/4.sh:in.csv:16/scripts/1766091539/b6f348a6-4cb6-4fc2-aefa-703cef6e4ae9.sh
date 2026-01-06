export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_vtZsjBm/ 
mkdir -p /tmp/pash_vtZsjBm/b6370c5edb40433ead435fd8384e5c43/ 
mkdir -p /tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/ 
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
{ rm -f "/tmp/pash_vtZsjBm/b6370c5edb40433ead435fd8384e5c43/#fifo96" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/b6370c5edb40433ead435fd8384e5c43/#fifo112" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo319" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo320" ; } 
 { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo381" ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_vtZsjBm/b6370c5edb40433ead435fd8384e5c43/#fifo96" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/b6370c5edb40433ead435fd8384e5c43/#fifo112" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo319" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo320" ; } 
 { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo381" ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_wrap bash -c ' cut -d "," -f 1 ' <"/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo320" >"/tmp/pash_vtZsjBm/b6370c5edb40433ead435fd8384e5c43/#fifo96" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_unwrap <"/tmp/pash_vtZsjBm/b6370c5edb40433ead435fd8384e5c43/#fifo96" >"/tmp/pash_vtZsjBm/b6370c5edb40433ead435fd8384e5c43/#fifo112" & }
pids_to_kill="${!} ${pids_to_kill}"
{ sort <"/tmp/pash_vtZsjBm/b6370c5edb40433ead435fd8384e5c43/#fifo112" >"/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo381" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo319" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo320" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*9ba391f4-52a6-419a-b1bd-95c2aae58586*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo319 send*d42ade91-e6c0-415e-bbf0-db7989833b3d*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo381 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
