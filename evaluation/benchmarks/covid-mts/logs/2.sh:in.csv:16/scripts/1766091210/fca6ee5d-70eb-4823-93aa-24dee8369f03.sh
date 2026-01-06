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
{ rm -f "/tmp/pash_cahH684/566f3e599dfa4a999b957742bec6928b/#fifo6" ; } 
 { { rm -f "/tmp/pash_cahH684/566f3e599dfa4a999b957742bec6928b/#fifo8" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo229" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo230" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo233" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo234" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo237" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo238" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo241" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo242" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo245" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo246" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo249" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo250" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo253" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo254" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo257" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo258" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo261" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo262" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo265" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo266" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo269" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo270" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo273" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo274" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo277" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo278" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo281" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo282" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo285" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo286" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo289" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo290" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo291" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo295" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo299" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo303" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo307" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo311" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo315" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo319" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo323" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo327" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo331" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo335" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo339" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo343" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo347" ; } 
 { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo351" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_cahH684/566f3e599dfa4a999b957742bec6928b/#fifo6" ; } 
 { { mkfifo "/tmp/pash_cahH684/566f3e599dfa4a999b957742bec6928b/#fifo8" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo229" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo230" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo233" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo234" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo237" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo238" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo241" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo242" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo245" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo246" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo249" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo250" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo253" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo254" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo257" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo258" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo261" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo262" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo265" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo266" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo269" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo270" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo273" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo274" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo277" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo278" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo281" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo282" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo285" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo286" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo289" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo290" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo291" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo295" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo299" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo303" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo307" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo311" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo315" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo319" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo323" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo327" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo331" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo335" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo339" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo343" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo347" ; } 
 { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo351" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_merge "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo230" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo234" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo238" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo242" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo246" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo250" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo254" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo258" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo262" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo266" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo270" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo274" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo278" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo282" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo286" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo290" >"/tmp/pash_cahH684/566f3e599dfa4a999b957742bec6928b/#fifo6" & }
pids_to_kill="${!} ${pids_to_kill}"
{ sort -u <"/tmp/pash_cahH684/566f3e599dfa4a999b957742bec6928b/#fifo6" >"/tmp/pash_cahH684/566f3e599dfa4a999b957742bec6928b/#fifo8" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_split "/tmp/pash_cahH684/566f3e599dfa4a999b957742bec6928b/#fifo8" 1000000 "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo291" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo295" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo299" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo303" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo307" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo311" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo315" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo319" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo323" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo327" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo331" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo335" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo339" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo343" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo347" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo351" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo229" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo230" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo233" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo234" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo237" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo238" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo241" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo242" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo245" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo246" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo249" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo250" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo253" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo254" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo257" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo258" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo261" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo262" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo265" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo266" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo269" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo270" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo273" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo274" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo277" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo278" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo281" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo282" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo285" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo286" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo289" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo290" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*537c7afc-2209-40fb-9c91-1169b46d1d2b*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo229 recv*bbbd9f71-dead-4938-821e-3b46ef103623*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo233 recv*4666bf20-0fbf-483e-8861-1ae4879e3671*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo237 recv*0a3fca0f-5fc4-4eda-a150-8306ff4eb623*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo241 recv*970ab3df-093a-46bb-a92c-ac263131eb43*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo245 recv*0758b0bb-1a73-49f7-bf2a-0e18836c7a20*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo249 recv*90e33c20-2a63-4005-9114-655c3d41b782*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo253 recv*962e27bf-cd0c-407b-934b-d04bf8ac016f*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo257 recv*0457e944-e281-4e98-a7a3-4ba955194360*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo261 recv*c5b20c1e-a79c-4c32-8fb0-a320636643dd*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo265 recv*3c793463-15ee-4d30-9b71-45f362794dc7*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo269 recv*3bbd87b5-28c1-4e0d-be72-43255f0801bd*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo273 recv*6328748d-1eb6-49ca-9ac3-595c05adfdec*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo277 recv*18e688f2-89df-4fe4-bdb3-2956b48c6e88*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo281 recv*591230f5-1aec-452f-b004-a133e9919ec6*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo285 recv*fea92475-5905-428f-a199-16af2f6b2c0b*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo289 send*8966b1ed-0f3c-4dcd-a3fa-a19dcb88b8b3*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo291 send*76c7648c-f239-41bb-96fd-8b993d4b4978*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo295 send*3ac90011-8c9a-4079-8758-282c2379183d*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo299 send*a405e37b-a274-41f9-8cb1-1399a1104831*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo303 send*07f59d49-364d-48ab-91f0-766739b1b184*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo307 send*f3472297-0d80-42a2-8953-7cdf7b650656*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo311 send*6031ddf0-4cba-440f-b9df-6c9d93044eb1*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo315 send*2973367e-bda3-48d6-8f47-5217756a81ab*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo319 send*0ca9a57f-b7ce-49b5-8f0e-9e4ddaf1cf50*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo323 send*4e803aba-f27b-4a9e-bffa-fdc5afa9255b*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo327 send*fd099f52-5f50-47c2-ac02-9686986b6d08*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo331 send*9a855c59-1579-4507-ab43-ed37f76aaed6*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo335 send*e6452a11-7d70-4349-89f8-1abbfa281082*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo339 send*2cf385fb-68c5-4e10-8688-050985f868c1*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo343 send*ebd3b993-307f-4a34-a9bf-516356bccb6c*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo347 send*ca3329f0-b201-4074-9983-4b6c3e0893b1*0*1*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo351 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
