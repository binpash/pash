export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_cahH684/ 
mkdir -p /tmp/pash_cahH684/566f3e599dfa4a999b957742bec6928b/ 
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
{ rm -f "/tmp/pash_cahH684/566f3e599dfa4a999b957742bec6928b/#fifo12" ; } 
 { { rm -f "/tmp/pash_cahH684/566f3e599dfa4a999b957742bec6928b/#fifo14" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo357" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo358" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo361" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo362" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo365" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo366" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo369" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo370" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo373" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo374" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo377" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo378" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo381" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo382" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo385" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo386" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo389" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo390" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo393" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo394" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo397" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo398" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo401" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo402" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo405" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo406" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo409" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo410" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo413" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo414" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo417" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo418" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo419" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo423" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo427" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo431" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo435" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo439" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo443" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo447" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo451" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo455" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo459" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo463" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo467" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo471" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo475" ; } 
 { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo479" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_cahH684/566f3e599dfa4a999b957742bec6928b/#fifo12" ; } 
 { { mkfifo "/tmp/pash_cahH684/566f3e599dfa4a999b957742bec6928b/#fifo14" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo357" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo358" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo361" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo362" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo365" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo366" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo369" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo370" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo373" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo374" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo377" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo378" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo381" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo382" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo385" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo386" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo389" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo390" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo393" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo394" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo397" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo398" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo401" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo402" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo405" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo406" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo409" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo410" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo413" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo414" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo417" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo418" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo419" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo423" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo427" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo431" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo435" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo439" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo443" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo447" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo451" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo455" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo459" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo463" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo467" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo471" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo475" ; } 
 { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo479" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ sort -m "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo358" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo362" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo366" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo370" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo374" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo378" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo382" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo386" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo390" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo394" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo398" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo402" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo406" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo410" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo414" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo418" >"/tmp/pash_cahH684/566f3e599dfa4a999b957742bec6928b/#fifo12" & }
pids_to_kill="${!} ${pids_to_kill}"
{ uniq -c <"/tmp/pash_cahH684/566f3e599dfa4a999b957742bec6928b/#fifo12" >"/tmp/pash_cahH684/566f3e599dfa4a999b957742bec6928b/#fifo14" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_split -r "/tmp/pash_cahH684/566f3e599dfa4a999b957742bec6928b/#fifo14" 1000000 "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo419" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo423" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo427" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo431" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo435" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo439" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo443" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo447" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo451" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo455" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo459" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo463" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo467" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo471" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo475" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo479" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo357" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo358" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo361" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo362" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo365" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo366" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo369" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo370" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo373" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo374" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo377" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo378" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo381" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo382" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo385" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo386" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo389" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo390" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo393" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo394" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo397" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo398" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo401" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo402" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo405" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo406" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo409" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo410" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo413" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo414" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo417" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo418" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*53433980-0ba2-48ae-8c2d-d074cbec7bdb*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo357 recv*1d9f3702-8f79-4f5a-b1f2-b4728a4f2c9d*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo361 recv*4785ad80-f410-465c-b8dd-d05b36452835*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo365 recv*553fdbdc-bb6a-4f71-8c5b-91698b72de9c*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo369 recv*9cc7c7db-78c5-4ec6-a8e3-4dee2f777f6d*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo373 recv*e7503297-20bd-4dbe-87c3-997bb99777d2*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo377 recv*0555a32d-842a-48fd-8373-3d1e761bf79d*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo381 recv*902a3a2f-20ec-4b0d-be55-a9caca838d88*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo385 recv*5d6d9dc2-c6e9-4982-a5a2-3e704fd345b9*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo389 recv*d4665b6c-858e-4e7b-b0a6-ef8d19bfc06b*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo393 recv*225bb67a-82d6-4f03-a57f-2b829177f780*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo397 recv*d08cacca-a8d8-4e74-9b30-8114c360fd20*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo401 recv*8f2cd59e-9aad-46b2-9137-09e706d69bd4*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo405 recv*592ac7ac-4592-4690-9854-7d8d3da86570*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo409 recv*0bf54d06-ef87-4200-8bb0-8c596240cab1*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo413 recv*076bdfef-aa8d-4607-913e-63cc506dd3a4*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo417 send*986f005a-79d1-40af-a146-f57d60a9db1d*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo419 send*7fc29827-d624-47f3-a356-d369524025a2*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo423 send*a3cfc986-66c6-4154-9f1a-66b97ca05ed2*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo427 send*9af4dae8-3f34-481e-ab00-f5dcc71c75a4*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo431 send*12dd9b08-2b69-44ac-891f-1e69c86272c8*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo435 send*b197eed6-af8d-4eb4-a95e-aec32677ec8e*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo439 send*09ca676b-49ab-40f9-9bdb-12cfd270d5ed*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo443 send*228495aa-9da7-4de0-ae50-301e50c98267*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo447 send*bd78b7eb-c9c5-48e2-a124-5e5006fbc4da*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo451 send*d22944ab-4869-4c8f-ae71-4fb7bc744df7*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo455 send*5b4969b0-fd60-48e4-ab30-e1c947baa22a*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo459 send*1fa82aac-3751-4d40-acf3-1714caff005d*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo463 send*1772c55b-7893-447e-ac80-6a666e715048*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo467 send*42f41751-9de9-447f-96ce-c9a2ca4e52e0*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo471 send*2a24a2b1-5415-41cf-ad68-1070e3795c38*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo475 send*35b527f7-febe-4c5e-b22a-5d169f8606ee*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo479 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
