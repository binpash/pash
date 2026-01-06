export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_gbcFrw6/ 
mkdir -p /tmp/pash_gbcFrw6/cd4ea6dd729a4293aa848eb873385dd9/ 
mkdir -p /tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/ 
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
{ rm -f "/tmp/pash_gbcFrw6/cd4ea6dd729a4293aa848eb873385dd9/#fifo6" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/cd4ea6dd729a4293aa848eb873385dd9/#fifo8" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo195" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo196" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo199" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo200" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo203" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo204" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo207" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo208" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo211" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo212" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo215" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo216" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo219" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo220" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo223" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo224" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo227" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo228" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo231" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo232" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo235" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo236" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo239" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo240" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo243" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo244" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo247" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo248" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo251" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo252" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo255" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo256" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo257" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo261" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo265" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo269" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo273" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo277" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo281" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo285" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo289" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo293" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo297" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo301" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo305" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo309" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo313" ; } 
 { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo317" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_gbcFrw6/cd4ea6dd729a4293aa848eb873385dd9/#fifo6" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/cd4ea6dd729a4293aa848eb873385dd9/#fifo8" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo195" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo196" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo199" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo200" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo203" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo204" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo207" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo208" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo211" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo212" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo215" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo216" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo219" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo220" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo223" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo224" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo227" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo228" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo231" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo232" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo235" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo236" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo239" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo240" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo243" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo244" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo247" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo248" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo251" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo252" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo255" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo256" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo257" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo261" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo265" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo269" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo273" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo277" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo281" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo285" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo289" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo293" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo297" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo301" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo305" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo309" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo313" ; } 
 { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo317" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_merge "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo196" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo200" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo204" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo208" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo212" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo216" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo220" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo224" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo228" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo232" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo236" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo240" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo244" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo248" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo252" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo256" >"/tmp/pash_gbcFrw6/cd4ea6dd729a4293aa848eb873385dd9/#fifo6" & }
pids_to_kill="${!} ${pids_to_kill}"
{ sort -u <"/tmp/pash_gbcFrw6/cd4ea6dd729a4293aa848eb873385dd9/#fifo6" >"/tmp/pash_gbcFrw6/cd4ea6dd729a4293aa848eb873385dd9/#fifo8" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_split "/tmp/pash_gbcFrw6/cd4ea6dd729a4293aa848eb873385dd9/#fifo8" 1000000 "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo257" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo261" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo265" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo269" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo273" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo277" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo281" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo285" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo289" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo293" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo297" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo301" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo305" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo309" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo313" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo317" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo195" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo196" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo199" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo200" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo203" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo204" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo207" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo208" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo211" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo212" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo215" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo216" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo219" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo220" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo223" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo224" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo227" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo228" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo231" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo232" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo235" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo236" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo239" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo240" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo243" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo244" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo247" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo248" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo251" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo252" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo255" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo256" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*47ac566c-805f-4a8a-bca1-54d02b3c613f*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo195 recv*1e8fac9f-5403-4256-9c5e-552d862c0cc3*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo199 recv*ad46b3b2-5c88-48e7-829e-9ea712ae9ae6*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo203 recv*6dd0889d-bad4-48d2-adcb-40ef7a3c0bb0*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo207 recv*9b19c60f-b6f1-4e51-a77b-103cd38bdda2*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo211 recv*73eae35f-9810-49fb-bbae-47bcfa37a0eb*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo215 recv*807fb9a0-f01d-4cac-9542-e1ddb39f882b*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo219 recv*ca9e8914-de6c-4d23-b4ea-c9f9ac86329b*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo223 recv*bc2e9968-efe0-4459-a6dc-466d9051b823*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo227 recv*19ad8e67-5306-47bd-bb7e-313ee37ec549*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo231 recv*bf96bbdc-cfec-4de2-b09b-3c13af27bed1*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo235 recv*9ce5fe55-cfbc-44b5-a8ac-0ff3bced8ec9*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo239 recv*5893dd79-ec66-4939-89f2-34d7fa25e72e*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo243 recv*5661052b-c1a6-4568-b840-de1205001613*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo247 recv*d84b8b10-bee1-4f85-9859-f6d5633151b1*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo251 recv*24d97713-2f87-48a6-b738-0727faf681cb*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo255 send*05f81314-d539-4976-81b0-c733b16a4888*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo257 send*2e8e4778-a037-4452-921e-c9ff9495d640*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo261 send*2c9c81b6-f436-456c-a753-6d823dc8bae7*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo265 send*7f98a473-b93d-43c4-b65a-fe0d54c4ceb0*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo269 send*28253fcf-811c-46c5-a3bd-d99ed5edbe24*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo273 send*135dad7c-7e7f-4edb-8a99-c2c0cbfb0b0f*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo277 send*3f543fc5-a431-446b-b414-6b8a8ec071f1*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo281 send*fa33ec19-a8a7-46f5-a09c-aa161a09368a*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo285 send*9cc43827-8163-4489-953f-040ebc7032b4*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo289 send*151dfa14-0c73-4bf4-8660-d412250fde3c*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo293 send*1ad9ab2b-10fd-4db4-b32f-ba3b4ea44a63*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo297 send*5c8881d7-285d-411e-8514-fb37f2346926*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo301 send*53430f99-b901-4ca0-8b15-f820eab9cef1*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo305 send*04484c08-60cc-4aaf-a36d-b6ccbbe15165*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo309 send*062ecfe1-6c9e-4339-9a57-71eaf3b0108b*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo313 send*603c4dca-45e2-4525-a697-a6dd518d742b*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo317 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
