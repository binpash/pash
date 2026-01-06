export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_psR0tRP/ 
mkdir -p /tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/ 
mkdir -p /tmp/pash_psR0tRP/8f01866ff5f743a8b65995a89637b3f5/ 
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
{ rm -f "/tmp/pash_psR0tRP/8f01866ff5f743a8b65995a89637b3f5/#fifo12" ; } 
 { { rm -f "/tmp/pash_psR0tRP/8f01866ff5f743a8b65995a89637b3f5/#fifo14" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo357" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo358" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo361" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo362" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo365" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo366" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo369" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo370" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo373" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo374" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo377" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo378" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo381" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo382" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo385" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo386" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo389" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo390" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo393" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo394" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo397" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo398" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo401" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo402" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo405" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo406" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo409" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo410" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo413" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo414" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo417" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo418" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo419" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo423" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo427" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo431" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo435" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo439" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo443" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo447" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo451" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo455" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo459" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo463" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo467" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo471" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo475" ; } 
 { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo479" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_psR0tRP/8f01866ff5f743a8b65995a89637b3f5/#fifo12" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/8f01866ff5f743a8b65995a89637b3f5/#fifo14" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo357" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo358" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo361" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo362" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo365" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo366" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo369" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo370" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo373" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo374" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo377" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo378" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo381" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo382" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo385" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo386" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo389" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo390" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo393" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo394" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo397" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo398" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo401" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo402" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo405" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo406" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo409" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo410" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo413" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo414" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo417" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo418" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo419" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo423" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo427" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo431" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo435" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo439" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo443" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo447" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo451" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo455" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo459" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo463" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo467" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo471" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo475" ; } 
 { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo479" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ sort -m "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo358" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo362" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo366" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo370" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo374" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo378" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo382" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo386" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo390" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo394" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo398" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo402" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo406" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo410" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo414" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo418" >"/tmp/pash_psR0tRP/8f01866ff5f743a8b65995a89637b3f5/#fifo12" & }
pids_to_kill="${!} ${pids_to_kill}"
{ uniq -c <"/tmp/pash_psR0tRP/8f01866ff5f743a8b65995a89637b3f5/#fifo12" >"/tmp/pash_psR0tRP/8f01866ff5f743a8b65995a89637b3f5/#fifo14" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_split -r "/tmp/pash_psR0tRP/8f01866ff5f743a8b65995a89637b3f5/#fifo14" 1000000 "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo419" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo423" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo427" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo431" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo435" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo439" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo443" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo447" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo451" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo455" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo459" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo463" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo467" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo471" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo475" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo479" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo357" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo358" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo361" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo362" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo365" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo366" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo369" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo370" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo373" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo374" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo377" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo378" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo381" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo382" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo385" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo386" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo389" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo390" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo393" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo394" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo397" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo398" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo401" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo402" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo405" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo406" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo409" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo410" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo413" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo414" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo417" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo418" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*b99cd80d-7cfc-48a0-b37f-1a507c0e68c3*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo357 recv*51c314d5-ae94-4b1e-a671-0c238ef130a9*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo361 recv*aab63d5d-1a28-4a22-9281-eb517982f5cd*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo365 recv*f7449c91-dcf5-48a2-a055-3f7b49b3ffa3*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo369 recv*9fce88a8-edbb-4d19-bf31-9ceaacf6c6a6*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo373 recv*f3e256fa-265e-44c4-a941-bc869057d407*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo377 recv*9e175d7c-119f-47e1-ba33-82634b6a2e10*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo381 recv*a3869464-a7a4-4a8d-a7bf-189484b06206*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo385 recv*a09a7622-c665-4986-b345-609bb0680325*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo389 recv*965d615b-68f0-42f0-9217-07d6ae98b5c6*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo393 recv*96aeca4a-8321-434d-97db-01bef9478e58*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo397 recv*d7875d2c-353d-4b33-80b6-1b0d6009dea4*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo401 recv*051ec091-ab18-46b7-8b21-ba833cf180e5*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo405 recv*7ba9a932-e19c-4c36-a07b-651594223d89*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo409 recv*c269e724-73dc-4504-b622-996ce28432e1*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo413 recv*2e8fe0c6-4c6d-4613-90e6-ed46661ddf89*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo417 send*48fd5897-9224-4bec-b5ba-fdac94b6f117*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo419 send*e46455be-f80d-49bd-9bd4-235b44f1686e*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo423 send*81501cd8-f893-48e3-ae9e-618bec1d51be*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo427 send*9b50bf83-0122-49c1-acc3-30de9973a6b9*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo431 send*6478fb1b-bc31-4ec5-a351-604dfab7ccbb*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo435 send*e774d543-c4bf-4436-bbbb-6aa5cd1f0b47*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo439 send*f53f75b9-9800-459c-8daf-e1b5dafd350c*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo443 send*5cecd62c-b090-44ac-8c00-1422e49246cb*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo447 send*9c8593c2-0b7a-4848-bff3-0d48e42e0cde*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo451 send*8bb67bc7-3629-4198-96ab-35de8ffbeced*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo455 send*67a3ff52-bf71-40c7-8b04-0d61771ac581*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo459 send*ea669c56-6574-4f99-aa2f-8b3882cc72c9*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo463 send*21565924-f86b-4ce1-916b-e2264a8ab028*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo467 send*6a2284c1-709c-40cf-a3cd-74718b8d1ae9*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo471 send*c823a332-63b1-4663-bf5b-c697a9d06cd9*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo475 send*2950a5b2-ca27-4a10-97b3-a7020ef835a2*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo479 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
