export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_vweH2q0/ 
mkdir -p /tmp/pash_vweH2q0/daa25b89ecd24ebaac69e20f741e9cc9/ 
mkdir -p /tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/ 
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
{ rm -f "/tmp/pash_vweH2q0/daa25b89ecd24ebaac69e20f741e9cc9/#fifo96" ; } 
 { { rm -f "/tmp/pash_vweH2q0/daa25b89ecd24ebaac69e20f741e9cc9/#fifo112" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo255" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo256" ; } 
 { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo317" ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_vweH2q0/daa25b89ecd24ebaac69e20f741e9cc9/#fifo96" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/daa25b89ecd24ebaac69e20f741e9cc9/#fifo112" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo255" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo256" ; } 
 { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo317" ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_wrap bash -c ' cut -d "," -f 1 ' <"/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo256" >"/tmp/pash_vweH2q0/daa25b89ecd24ebaac69e20f741e9cc9/#fifo96" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_unwrap <"/tmp/pash_vweH2q0/daa25b89ecd24ebaac69e20f741e9cc9/#fifo96" >"/tmp/pash_vweH2q0/daa25b89ecd24ebaac69e20f741e9cc9/#fifo112" & }
pids_to_kill="${!} ${pids_to_kill}"
{ sort <"/tmp/pash_vweH2q0/daa25b89ecd24ebaac69e20f741e9cc9/#fifo112" >"/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo317" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo255" -o "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo256" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*a1fa1f26-3ad6-4ec3-8d49-f03e95f2369a*1*0*/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo255 send*ac79a8e3-2a86-4765-bfa2-c04df23b4acd*0*1*/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo317 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
