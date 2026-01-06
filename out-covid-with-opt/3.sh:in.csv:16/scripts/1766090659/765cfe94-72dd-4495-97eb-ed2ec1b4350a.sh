export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_oJI9zDr/ 
mkdir -p /tmp/pash_oJI9zDr/ceecddd825db445bab2b7187bbc92eb0/ 
mkdir -p /tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/ 
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
{ rm -f "/tmp/pash_oJI9zDr/ceecddd825db445bab2b7187bbc92eb0/#fifo6" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/ceecddd825db445bab2b7187bbc92eb0/#fifo8" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo165" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo166" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo169" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo170" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo173" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo174" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo177" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo178" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo181" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo182" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo185" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo186" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo189" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo190" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo193" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo194" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo197" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo198" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo201" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo202" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo205" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo206" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo209" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo210" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo213" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo214" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo217" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo218" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo221" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo222" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo225" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo226" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo227" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo231" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo235" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo239" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo243" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo247" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo251" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo255" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo259" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo263" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo267" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo271" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo275" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo279" ; } 
 { { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo283" ; } 
 { rm -f "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo287" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_oJI9zDr/ceecddd825db445bab2b7187bbc92eb0/#fifo6" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/ceecddd825db445bab2b7187bbc92eb0/#fifo8" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo165" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo166" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo169" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo170" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo173" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo174" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo177" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo178" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo181" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo182" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo185" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo186" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo189" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo190" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo193" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo194" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo197" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo198" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo201" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo202" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo205" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo206" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo209" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo210" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo213" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo214" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo217" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo218" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo221" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo222" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo225" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo226" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo227" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo231" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo235" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo239" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo243" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo247" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo251" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo255" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo259" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo263" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo267" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo271" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo275" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo279" ; } 
 { { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo283" ; } 
 { mkfifo "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo287" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_merge "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo166" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo170" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo174" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo178" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo182" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo186" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo190" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo194" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo198" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo202" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo206" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo210" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo214" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo218" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo222" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo226" >"/tmp/pash_oJI9zDr/ceecddd825db445bab2b7187bbc92eb0/#fifo6" & }
pids_to_kill="${!} ${pids_to_kill}"
{ sort -u <"/tmp/pash_oJI9zDr/ceecddd825db445bab2b7187bbc92eb0/#fifo6" >"/tmp/pash_oJI9zDr/ceecddd825db445bab2b7187bbc92eb0/#fifo8" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_split "/tmp/pash_oJI9zDr/ceecddd825db445bab2b7187bbc92eb0/#fifo8" 1000000 "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo227" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo231" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo235" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo239" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo243" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo247" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo251" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo255" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo259" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo263" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo267" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo271" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo275" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo279" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo283" "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo287" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo165" -o "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo166" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo169" -o "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo170" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo173" -o "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo174" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo177" -o "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo178" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo181" -o "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo182" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo185" -o "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo186" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo189" -o "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo190" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo193" -o "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo194" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo197" -o "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo198" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo201" -o "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo202" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo205" -o "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo206" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo209" -o "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo210" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo213" -o "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo214" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo217" -o "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo218" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo221" -o "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo222" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo225" -o "/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo226" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*cbd099c0-ff32-4f2d-b69d-dd1617447536*1*0*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo165 recv*b546884a-d1da-43eb-b0ec-7c0b07d2f0e7*1*0*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo169 recv*bbfc07d5-cdf3-4095-a12b-582317b40627*1*0*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo173 recv*f7a58e8e-b241-4654-b17f-5910f20210d5*1*0*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo177 recv*5c926af6-d52c-46d1-aeab-2e44ff788ba6*1*0*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo181 recv*a69d3c9f-2b92-460a-8107-e5371ceffe8b*1*0*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo185 recv*e91008e8-a31b-4e3e-9705-f5995da486ca*1*0*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo189 recv*19a6db65-2297-47c1-847c-f9befbbb0fd2*1*0*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo193 recv*57c2c3c7-62dc-4716-bfc2-1db26f40e7e1*1*0*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo197 recv*1f1020a8-69c5-487c-8bc0-3a8a98358e08*1*0*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo201 recv*3ee4fc52-e50a-4ad7-9ca2-bcff00f1eb46*1*0*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo205 recv*e15600db-7e1e-4ee0-b9a8-62c9b8284952*1*0*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo209 recv*a9fa9348-7d6c-4d08-96ec-3a944eb53bff*1*0*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo213 recv*d8501416-693f-4f18-8a4a-db0f45fccfdc*1*0*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo217 recv*fe21c310-674f-496c-88dd-7218f123dbba*1*0*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo221 recv*d1b1f938-bbc8-4f7b-8cfd-268dc480d55a*1*0*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo225 send*931029b9-348b-4065-9fc6-a272ead655c1*0*1*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo227 send*0da6d4e1-dbdf-4dfd-addb-f5217ca7a7f1*0*1*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo231 send*db373801-c331-48aa-a8d9-c2906a1c466d*0*1*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo235 send*f2acbd34-7ffd-4adf-acb1-27d9bc208d1e*0*1*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo239 send*7d61757f-dd06-49d1-a174-bcd1c78422cd*0*1*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo243 send*1c3b9e07-cb3b-4402-8796-b65721cac896*0*1*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo247 send*8c2cdc8f-0295-4e9f-98e6-d5081fbe1e61*0*1*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo251 send*de0cd902-9f30-4769-9d66-26309ec54619*0*1*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo255 send*f0e15554-e89e-4066-a0b3-929f357f6d8b*0*1*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo259 send*6a6ddf3a-b143-4386-8c54-ada0469fe041*0*1*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo263 send*15165b6c-6614-4f25-b5dc-7bbdf0fb71c4*0*1*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo267 send*54b9c773-9a86-42cd-a7de-0c1136b9cb0d*0*1*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo271 send*e246ff32-2d00-4b8b-9573-8162059197ff*0*1*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo275 send*9cb7d59d-16c3-4149-afd1-8167dcfd4e01*0*1*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo279 send*f8e94aab-3e8d-47f3-b190-9da302f324ae*0*1*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo283 send*2ca490e4-896f-49d9-8a73-d7384fc92ae5*0*1*/tmp/pash_oJI9zDr/fa2afb37cd0543cbb5628bd247b9b3ae/#fifo287 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
