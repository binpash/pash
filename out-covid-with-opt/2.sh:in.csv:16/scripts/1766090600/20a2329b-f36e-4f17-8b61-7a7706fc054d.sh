export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_k1vQSpy/ 
mkdir -p /tmp/pash_k1vQSpy/32e98ac9d8d643d6ac226bf864b2a19e/ 
mkdir -p /tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/ 
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
{ rm -f "/tmp/pash_k1vQSpy/32e98ac9d8d643d6ac226bf864b2a19e/#fifo6" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/32e98ac9d8d643d6ac226bf864b2a19e/#fifo8" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo165" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo166" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo169" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo170" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo173" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo174" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo177" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo178" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo181" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo182" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo185" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo186" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo189" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo190" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo193" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo194" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo197" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo198" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo201" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo202" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo205" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo206" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo209" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo210" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo213" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo214" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo217" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo218" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo221" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo222" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo225" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo226" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo227" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo231" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo235" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo239" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo243" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo247" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo251" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo255" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo259" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo263" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo267" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo271" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo275" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo279" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo283" ; } 
 { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo287" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_k1vQSpy/32e98ac9d8d643d6ac226bf864b2a19e/#fifo6" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/32e98ac9d8d643d6ac226bf864b2a19e/#fifo8" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo165" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo166" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo169" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo170" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo173" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo174" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo177" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo178" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo181" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo182" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo185" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo186" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo189" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo190" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo193" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo194" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo197" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo198" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo201" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo202" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo205" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo206" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo209" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo210" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo213" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo214" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo217" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo218" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo221" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo222" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo225" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo226" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo227" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo231" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo235" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo239" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo243" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo247" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo251" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo255" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo259" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo263" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo267" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo271" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo275" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo279" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo283" ; } 
 { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo287" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_merge "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo166" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo170" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo174" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo178" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo182" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo186" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo190" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo194" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo198" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo202" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo206" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo210" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo214" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo218" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo222" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo226" >"/tmp/pash_k1vQSpy/32e98ac9d8d643d6ac226bf864b2a19e/#fifo6" & }
pids_to_kill="${!} ${pids_to_kill}"
{ sort -u <"/tmp/pash_k1vQSpy/32e98ac9d8d643d6ac226bf864b2a19e/#fifo6" >"/tmp/pash_k1vQSpy/32e98ac9d8d643d6ac226bf864b2a19e/#fifo8" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_split "/tmp/pash_k1vQSpy/32e98ac9d8d643d6ac226bf864b2a19e/#fifo8" 1000000 "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo227" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo231" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo235" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo239" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo243" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo247" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo251" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo255" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo259" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo263" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo267" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo271" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo275" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo279" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo283" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo287" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo165" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo166" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo169" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo170" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo173" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo174" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo177" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo178" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo181" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo182" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo185" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo186" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo189" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo190" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo193" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo194" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo197" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo198" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo201" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo202" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo205" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo206" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo209" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo210" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo213" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo214" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo217" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo218" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo221" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo222" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo225" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo226" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*5143e30d-43ac-451a-88ff-91b70a6764bb*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo165 recv*e492d916-c4c5-46a8-bbfe-e81785a9af66*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo169 recv*5b5b90b0-e6f8-41fb-b940-1e6b4de7624f*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo173 recv*ae5612a6-ec06-4574-9835-ecb0b6b83739*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo177 recv*9bb51b82-4cb7-47d9-afe3-28c98d8ef9a3*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo181 recv*e09ee467-58dc-4a22-bc67-4f654997e238*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo185 recv*ca13fac0-66d4-4265-ad9d-8606cfc55ffa*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo189 recv*539f75af-b0dd-4131-afd0-f2d58c0d8bd9*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo193 recv*d1b38654-b31a-4282-8687-f81e941b8eb0*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo197 recv*a6f8de7f-6386-41e2-a439-85f6c3167f96*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo201 recv*d5506763-7df8-4aff-9205-ebe807d6126c*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo205 recv*96ccc90b-56a7-490b-8f91-255a00bedec9*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo209 recv*2a65123c-5bbd-49e2-881c-719c1283d307*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo213 recv*17552737-cd20-4811-aa69-bb1aff7eec11*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo217 recv*98e5d71a-ae3a-45a3-8a0a-8cee711ac4fe*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo221 recv*f3c5ff1a-7966-4014-878f-e2edfad74dfe*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo225 send*a75ae3cc-0ae4-4fc6-93d9-430fc6cba50b*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo227 send*94e00054-fbe3-4f24-9c1c-b4c8a65d63ad*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo231 send*512b734c-6f99-471d-953b-94b671a9a70b*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo235 send*2f797561-e06c-4439-9b95-75f41d5b645a*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo239 send*6a16593b-9964-4620-91b2-8961a468a1ff*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo243 send*0d6eafe0-3839-406a-93c2-b5606876da99*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo247 send*3caf595b-cbe7-4244-ada6-60dcb0d84476*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo251 send*3ee73516-f2a9-447d-af04-7b97aa6d3e52*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo255 send*fd8dedc3-4673-4ccf-bab5-0217cc2db440*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo259 send*2490aa65-2583-443d-8e0f-f25c4534b316*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo263 send*e53d17db-2589-422d-a4c2-193e0c865c80*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo267 send*e8ce0819-e147-49de-af3c-996c07861ff7*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo271 send*78601fda-e869-4db0-a410-bb1e3aa05356*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo275 send*280c10bd-3aec-4af1-9f44-266beb2bb679*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo279 send*f8f8fbba-ec06-40ca-beaa-cad35fa66f8f*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo283 send*f68ce6ca-8a9f-40cb-b721-a827d0e6c6e1*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo287 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
