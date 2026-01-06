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
{ rm -f "/tmp/pash_psR0tRP/8f01866ff5f743a8b65995a89637b3f5/#fifo6" ; } 
 { { rm -f "/tmp/pash_psR0tRP/8f01866ff5f743a8b65995a89637b3f5/#fifo8" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo229" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo230" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo233" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo234" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo237" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo238" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo241" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo242" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo245" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo246" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo249" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo250" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo253" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo254" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo257" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo258" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo261" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo262" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo265" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo266" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo269" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo270" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo273" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo274" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo277" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo278" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo281" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo282" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo285" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo286" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo289" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo290" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo291" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo295" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo299" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo303" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo307" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo311" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo315" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo319" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo323" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo327" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo331" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo335" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo339" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo343" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo347" ; } 
 { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo351" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_psR0tRP/8f01866ff5f743a8b65995a89637b3f5/#fifo6" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/8f01866ff5f743a8b65995a89637b3f5/#fifo8" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo229" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo230" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo233" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo234" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo237" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo238" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo241" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo242" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo245" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo246" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo249" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo250" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo253" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo254" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo257" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo258" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo261" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo262" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo265" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo266" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo269" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo270" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo273" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo274" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo277" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo278" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo281" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo282" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo285" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo286" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo289" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo290" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo291" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo295" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo299" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo303" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo307" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo311" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo315" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo319" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo323" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo327" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo331" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo335" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo339" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo343" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo347" ; } 
 { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo351" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_merge "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo230" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo234" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo238" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo242" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo246" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo250" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo254" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo258" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo262" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo266" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo270" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo274" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo278" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo282" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo286" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo290" >"/tmp/pash_psR0tRP/8f01866ff5f743a8b65995a89637b3f5/#fifo6" & }
pids_to_kill="${!} ${pids_to_kill}"
{ sort -u <"/tmp/pash_psR0tRP/8f01866ff5f743a8b65995a89637b3f5/#fifo6" >"/tmp/pash_psR0tRP/8f01866ff5f743a8b65995a89637b3f5/#fifo8" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_split "/tmp/pash_psR0tRP/8f01866ff5f743a8b65995a89637b3f5/#fifo8" 1000000 "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo291" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo295" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo299" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo303" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo307" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo311" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo315" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo319" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo323" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo327" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo331" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo335" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo339" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo343" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo347" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo351" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo229" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo230" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo233" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo234" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo237" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo238" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo241" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo242" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo245" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo246" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo249" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo250" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo253" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo254" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo257" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo258" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo261" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo262" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo265" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo266" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo269" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo270" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo273" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo274" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo277" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo278" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo281" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo282" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo285" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo286" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo289" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo290" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*791462da-1405-4295-91cf-63d1417a41bc*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo229 recv*109e2624-44b0-4811-800c-cb996bde6bb6*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo233 recv*765ddb73-53e8-4cf1-96e3-6c0698466b4d*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo237 recv*eb0805f2-0fb3-48f7-a584-f476af8095c7*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo241 recv*05aa12f0-a496-4819-9048-2bbb2ff3c341*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo245 recv*e8fa7897-feb8-4005-8a4a-57ded254fc46*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo249 recv*145ae13a-3a88-49f0-b992-e253f4563725*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo253 recv*0e9baf95-902e-46ec-865e-82ae112e826c*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo257 recv*8c12fb6d-7e02-4b9c-befa-33cd46d17c8c*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo261 recv*faf9143d-8f7b-4a3d-b468-a29f385f80f3*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo265 recv*a56c4a6a-aa96-4ec4-a655-92eda0f881ab*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo269 recv*0e64d3e2-6c94-4d64-93c3-f1272ef63994*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo273 recv*0a3a9c3c-8197-40b7-b6c2-bef22ca28ae3*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo277 recv*da2d7d90-21b7-4918-978c-b031f73be163*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo281 recv*74bc5ea4-ec73-40fa-9656-834c953cae6d*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo285 recv*80b16661-de3f-4d41-9fa4-aa34502fd412*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo289 send*fc9f9890-cc89-4b64-bdc4-b5b8bda6413d*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo291 send*bb7ba015-dbca-496f-b854-cfe2f6426fa2*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo295 send*bd667bdd-b423-4ab0-a3e2-b2da494cd2b7*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo299 send*969cc029-ae77-4a69-a97d-2c581e6f4a33*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo303 send*b7eff960-47dc-43e6-a3b7-81860ca89000*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo307 send*9769a005-d679-464b-9bc4-9a0e43a640c4*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo311 send*2aeb0278-c7f7-4508-bce7-cc6045137f20*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo315 send*682761c7-0dda-4f39-8af3-f6cc7c79f4e3*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo319 send*a028a240-c255-40c5-b10a-68376ba801e1*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo323 send*203ea626-6c42-4eb2-8c44-2dd6f8606a7f*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo327 send*ebccbccb-42fe-4326-a64e-cb893888d556*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo331 send*3a615116-ca46-4869-bd1e-583af6f1abf6*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo335 send*6bf28841-a497-46c1-ba77-812a923bb136*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo339 send*4a7e3c16-85b4-4875-8221-5301567840d2*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo343 send*506f4169-515f-41b5-99ac-69afacdaec86*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo347 send*6fd869d1-0c78-4fda-b1a8-2c9ed640fc7d*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo351 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
