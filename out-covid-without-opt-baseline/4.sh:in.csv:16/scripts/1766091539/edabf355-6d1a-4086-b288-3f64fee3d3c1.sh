export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_vtZsjBm/ 
mkdir -p /tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/ 
mkdir -p /tmp/pash_vtZsjBm/68a7f4eeec164679bf1bc1313be933aa/ 
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
{ rm -f "/tmp/pash_vtZsjBm/68a7f4eeec164679bf1bc1313be933aa/#fifo6" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/68a7f4eeec164679bf1bc1313be933aa/#fifo8" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo195" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo196" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo199" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo200" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo203" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo204" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo207" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo208" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo211" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo212" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo215" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo216" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo219" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo220" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo223" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo224" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo227" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo228" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo231" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo232" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo235" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo236" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo239" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo240" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo243" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo244" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo247" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo248" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo251" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo252" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo255" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo256" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo257" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo261" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo265" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo269" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo273" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo277" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo281" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo285" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo289" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo293" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo297" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo301" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo305" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo309" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo313" ; } 
 { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo317" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_vtZsjBm/68a7f4eeec164679bf1bc1313be933aa/#fifo6" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/68a7f4eeec164679bf1bc1313be933aa/#fifo8" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo195" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo196" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo199" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo200" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo203" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo204" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo207" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo208" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo211" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo212" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo215" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo216" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo219" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo220" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo223" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo224" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo227" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo228" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo231" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo232" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo235" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo236" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo239" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo240" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo243" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo244" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo247" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo248" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo251" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo252" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo255" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo256" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo257" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo261" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo265" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo269" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo273" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo277" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo281" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo285" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo289" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo293" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo297" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo301" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo305" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo309" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo313" ; } 
 { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo317" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_merge "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo196" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo200" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo204" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo208" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo212" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo216" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo220" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo224" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo228" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo232" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo236" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo240" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo244" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo248" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo252" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo256" >"/tmp/pash_vtZsjBm/68a7f4eeec164679bf1bc1313be933aa/#fifo6" & }
pids_to_kill="${!} ${pids_to_kill}"
{ sort -u <"/tmp/pash_vtZsjBm/68a7f4eeec164679bf1bc1313be933aa/#fifo6" >"/tmp/pash_vtZsjBm/68a7f4eeec164679bf1bc1313be933aa/#fifo8" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_split "/tmp/pash_vtZsjBm/68a7f4eeec164679bf1bc1313be933aa/#fifo8" 1000000 "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo257" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo261" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo265" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo269" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo273" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo277" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo281" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo285" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo289" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo293" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo297" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo301" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo305" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo309" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo313" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo317" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo195" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo196" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo199" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo200" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo203" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo204" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo207" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo208" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo211" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo212" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo215" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo216" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo219" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo220" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo223" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo224" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo227" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo228" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo231" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo232" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo235" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo236" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo239" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo240" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo243" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo244" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo247" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo248" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo251" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo252" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo255" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo256" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*b1da0abd-af21-4afd-922c-48d64e3fc188*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo195 recv*a1ef1839-6963-438b-bdec-862bce44d2c3*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo199 recv*98807c3e-ecb2-4165-83e7-551ef5d56dd4*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo203 recv*299da538-5187-4e74-b481-5edd37690536*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo207 recv*f155dff0-63ae-40ff-80d2-56712b7eb8e2*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo211 recv*ec4906d6-5904-43e8-b37e-c4636b583bd4*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo215 recv*25a6add4-b44c-40f7-95f0-06e6a9246e0c*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo219 recv*9e4a3de5-173c-4c9b-9e93-558d0b845549*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo223 recv*1da63ef1-94a7-4fc1-87ed-20be673a0bb4*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo227 recv*d8583e13-ca52-4404-915d-eff486e898e0*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo231 recv*7cccbcb3-b8c0-4dd2-b87c-9626c8238a2e*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo235 recv*720ec09d-b892-4bcf-af17-0ca606873564*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo239 recv*a3b7f3b1-5e30-41e4-96f4-8ef30187a4c0*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo243 recv*cfa83162-7f64-48d2-818c-87680b2bb845*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo247 recv*ff5cb20f-0977-4eeb-a01d-87292a074a16*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo251 recv*305cf441-707d-43b1-ba27-bdf11b5bbe84*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo255 send*446c6296-1e86-4888-89b2-636d68b6e616*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo257 send*f5e8c487-5533-47fc-8c5a-b96f35af10cc*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo261 send*24f0d5bd-aaae-42c1-a42c-85bfb37f491a*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo265 send*396ea422-94d8-4eb9-8402-7cacdb020d05*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo269 send*76cb3d8f-45c7-474d-b0ab-99a4a56377c3*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo273 send*2b2a1e10-a33b-41a6-a55f-96a4d4c9f056*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo277 send*4ff39bd8-784a-437a-9805-3e35a94a1854*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo281 send*bab4f15c-fe24-4b1d-86a3-220564a664f9*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo285 send*7634f722-709c-4b3d-8847-05cf8c9c88a2*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo289 send*2047ea10-372d-41dc-9a31-22f1c154a021*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo293 send*9baa4894-0e31-49bd-bb0d-138e560c5a99*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo297 send*c7d8c194-2577-4e51-ac15-4fdce7c5deed*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo301 send*20bea023-7cac-4dd1-8310-3539390c630c*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo305 send*b9f32ae9-5dde-4df4-9246-6791fa973d86*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo309 send*78ca56d2-a8b3-44f5-abe7-e88c9afa9561*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo313 send*9ba391f4-52a6-419a-b1bd-95c2aae58586*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo317 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
