export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_qrxfus7/ 
mkdir -p /tmp/pash_qrxfus7/9f66d8e02d6a4813aed714a670e9f9a1/ 
mkdir -p /tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/ 
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
{ rm -f "/tmp/pash_qrxfus7/9f66d8e02d6a4813aed714a670e9f9a1/#fifo6" ; } 
 { { rm -f "/tmp/pash_qrxfus7/9f66d8e02d6a4813aed714a670e9f9a1/#fifo8" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo131" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo132" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo135" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo136" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo139" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo140" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo143" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo144" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo147" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo148" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo151" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo152" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo155" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo156" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo159" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo160" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo163" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo164" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo167" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo168" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo171" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo172" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo175" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo176" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo179" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo180" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo183" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo184" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo187" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo188" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo191" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo192" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo193" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo197" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo201" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo205" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo209" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo213" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo217" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo221" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo225" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo229" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo233" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo237" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo241" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo245" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo249" ; } 
 { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo253" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_qrxfus7/9f66d8e02d6a4813aed714a670e9f9a1/#fifo6" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/9f66d8e02d6a4813aed714a670e9f9a1/#fifo8" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo131" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo132" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo135" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo136" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo139" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo140" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo143" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo144" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo147" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo148" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo151" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo152" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo155" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo156" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo159" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo160" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo163" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo164" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo167" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo168" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo171" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo172" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo175" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo176" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo179" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo180" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo183" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo184" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo187" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo188" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo191" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo192" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo193" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo197" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo201" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo205" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo209" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo213" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo217" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo221" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo225" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo229" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo233" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo237" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo241" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo245" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo249" ; } 
 { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo253" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_merge "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo132" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo136" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo140" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo144" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo148" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo152" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo156" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo160" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo164" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo168" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo172" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo176" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo180" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo184" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo188" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo192" >"/tmp/pash_qrxfus7/9f66d8e02d6a4813aed714a670e9f9a1/#fifo6" & }
pids_to_kill="${!} ${pids_to_kill}"
{ sort -u <"/tmp/pash_qrxfus7/9f66d8e02d6a4813aed714a670e9f9a1/#fifo6" >"/tmp/pash_qrxfus7/9f66d8e02d6a4813aed714a670e9f9a1/#fifo8" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_split "/tmp/pash_qrxfus7/9f66d8e02d6a4813aed714a670e9f9a1/#fifo8" 1000000 "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo193" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo197" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo201" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo205" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo209" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo213" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo217" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo221" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo225" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo229" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo233" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo237" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo241" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo245" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo249" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo253" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo131" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo132" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo135" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo136" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo139" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo140" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo143" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo144" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo147" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo148" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo151" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo152" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo155" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo156" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo159" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo160" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo163" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo164" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo167" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo168" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo171" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo172" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo175" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo176" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo179" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo180" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo183" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo184" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo187" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo188" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo191" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo192" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*eb62e4e4-2cb0-4e81-ba1d-717a2b8c355d*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo131 recv*d41c5959-6d8b-4e77-9834-b939180d8cb0*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo135 recv*93e44592-7c4d-44c1-87c9-4c5b87991356*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo139 recv*b123349b-5de4-426b-8be5-6c07d8c4056e*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo143 recv*2cf910b9-2c67-42af-ae66-e07529d566a2*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo147 recv*af7f8069-0d2e-450f-bb57-f4e93b4d4a08*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo151 recv*b636bf3c-6b07-4d62-a6a0-2c33ca661b7b*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo155 recv*9832e5aa-07da-4c39-bef4-ac7ab06460e0*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo159 recv*c7150e75-9327-48c3-aec2-9f30abc42551*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo163 recv*c0116ad9-3b34-4373-a9fa-b93f1f1f22e6*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo167 recv*7c987ee8-b54c-48cd-b765-938f974257a8*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo171 recv*fff49835-d2a8-4eb6-b7de-f60fed162284*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo175 recv*99e5c6e4-e5ce-4625-8902-2e4a5c20db96*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo179 recv*7376904c-4f23-42a1-85cd-a3e06b9060ee*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo183 recv*8ef91456-0be9-41ad-a591-8890b7c1fe6e*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo187 recv*059da216-fa42-429c-83b5-d1ec531ba189*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo191 send*00820ebe-36d9-46f3-87dd-d2396b309e27*0*1*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo193 send*6f3b27d3-067a-4509-8dfc-754be09b5dc2*0*1*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo197 send*30a80932-2751-4611-ba65-a070ccfb71c6*0*1*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo201 send*5c537cd6-a7f9-422b-87af-42405aaa56d2*0*1*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo205 send*af02c1e9-6148-4cd4-b249-336fde62f4fd*0*1*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo209 send*81b9bd45-80e5-4662-a6e6-15ceebe40bf3*0*1*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo213 send*89cb8ae5-8282-4fbd-95a0-022b10489839*0*1*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo217 send*270a54dc-57d8-40f2-aa29-42d3951e9c5d*0*1*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo221 send*313099e5-6e5b-4217-8a84-92177e5837b1*0*1*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo225 send*676895ce-d2b5-426c-925b-17dd9a70a2bc*0*1*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo229 send*82f9b09a-7ab4-4bc7-9aa8-51f8c1ba47cd*0*1*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo233 send*4a3eee8f-9d5e-46a1-bb5f-96226b3d5f8c*0*1*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo237 send*62d6f117-f6d3-4816-8e6c-3ba6a44b11c8*0*1*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo241 send*16f7d37a-93a5-413e-9591-720684358b03*0*1*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo245 send*dac75ef7-8399-473d-a987-567a4bec2061*0*1*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo249 send*73adf0cd-250b-4010-9992-54f40a48d88b*0*1*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo253 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
