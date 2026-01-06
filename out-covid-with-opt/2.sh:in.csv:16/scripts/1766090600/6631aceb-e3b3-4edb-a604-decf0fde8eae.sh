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
{ rm -f "/tmp/pash_k1vQSpy/32e98ac9d8d643d6ac226bf864b2a19e/#fifo12" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/32e98ac9d8d643d6ac226bf864b2a19e/#fifo14" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo293" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo294" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo297" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo298" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo301" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo302" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo305" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo306" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo309" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo310" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo313" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo314" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo317" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo318" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo321" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo322" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo325" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo326" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo329" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo330" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo333" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo334" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo337" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo338" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo341" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo342" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo345" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo346" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo349" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo350" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo353" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo354" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo355" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo359" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo363" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo367" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo371" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo375" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo379" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo383" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo387" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo391" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo395" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo399" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo403" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo407" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo411" ; } 
 { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo415" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_k1vQSpy/32e98ac9d8d643d6ac226bf864b2a19e/#fifo12" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/32e98ac9d8d643d6ac226bf864b2a19e/#fifo14" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo293" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo294" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo297" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo298" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo301" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo302" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo305" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo306" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo309" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo310" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo313" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo314" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo317" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo318" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo321" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo322" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo325" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo326" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo329" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo330" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo333" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo334" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo337" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo338" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo341" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo342" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo345" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo346" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo349" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo350" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo353" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo354" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo355" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo359" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo363" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo367" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo371" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo375" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo379" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo383" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo387" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo391" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo395" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo399" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo403" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo407" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo411" ; } 
 { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo415" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ sort -m "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo294" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo298" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo302" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo306" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo310" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo314" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo318" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo322" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo326" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo330" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo334" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo338" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo342" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo346" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo350" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo354" >"/tmp/pash_k1vQSpy/32e98ac9d8d643d6ac226bf864b2a19e/#fifo12" & }
pids_to_kill="${!} ${pids_to_kill}"
{ uniq -c <"/tmp/pash_k1vQSpy/32e98ac9d8d643d6ac226bf864b2a19e/#fifo12" >"/tmp/pash_k1vQSpy/32e98ac9d8d643d6ac226bf864b2a19e/#fifo14" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_split -r "/tmp/pash_k1vQSpy/32e98ac9d8d643d6ac226bf864b2a19e/#fifo14" 1000000 "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo355" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo359" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo363" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo367" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo371" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo375" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo379" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo383" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo387" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo391" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo395" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo399" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo403" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo407" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo411" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo415" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo293" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo294" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo297" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo298" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo301" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo302" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo305" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo306" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo309" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo310" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo313" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo314" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo317" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo318" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo321" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo322" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo325" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo326" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo329" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo330" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo333" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo334" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo337" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo338" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo341" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo342" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo345" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo346" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo349" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo350" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo353" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo354" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*bd386264-53e3-4e44-b3fc-1c9c600ba766*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo293 recv*c03dc005-d459-4300-9eb5-28a35d5007b0*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo297 recv*1f5541db-4d37-404d-a5be-341c65dba2db*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo301 recv*a7f71d62-5dfb-450e-b9e9-e6ac6598071c*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo305 recv*a1246dd7-dfab-4f69-8e69-383ee63b6897*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo309 recv*228a74c8-a596-4e81-b84d-3ee0b9076b7d*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo313 recv*1d9d0dcd-e761-428e-910b-34b9ff0630c6*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo317 recv*5a797090-72e6-49ae-9dbc-e4b137e0a5be*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo321 recv*8e1bc118-ddfe-4470-a5a9-a344c310ba26*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo325 recv*de5cfc89-5c9d-4b6f-8407-1f0b031de7df*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo329 recv*15395857-523e-45c7-baeb-98c24b45a662*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo333 recv*4a775f3f-79f5-486a-b459-b8f8f33afee6*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo337 recv*1678cd8b-a65f-4035-a260-025050c4ba75*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo341 recv*4db4b9a9-1ebf-4b94-bc53-74938a77d6cc*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo345 recv*af1ceedc-f516-4c98-be39-d20107b450e2*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo349 recv*27a56c89-d241-453a-b4b8-6fd2d85272b7*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo353 send*2da95de7-6372-4632-a68a-7cb22e7a99d1*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo355 send*8c20dd05-9f3e-4403-b43a-e27b9738dfc2*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo359 send*a50a3e1d-ea42-4f58-8bdb-cdc62830954e*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo363 send*33a90266-329f-4cd2-ab03-970e70c00ae9*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo367 send*2a8e7ecc-3328-4cab-85f1-f3ce2685a249*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo371 send*3632fb6d-a8c0-4988-a42f-0d3d14f19916*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo375 send*26fb9038-e2d1-43e0-a4eb-dfd3aee65526*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo379 send*c967e9dd-2f54-484d-b3cf-0330edf9feaf*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo383 send*3db06f16-427e-47fa-bfe2-a72b25c12bb5*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo387 send*c3e0dbe0-ba8f-4cb5-9c1f-de40ea5a4d9d*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo391 send*419bd2e1-083b-4df7-a702-f25b79f15ec6*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo395 send*dcea598a-cdf2-424a-b09f-c1d74161183a*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo399 send*a305d2a1-3dca-42f5-8bd6-a68ad8f69284*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo403 send*cd0df912-50a1-41ee-904c-046774648ec7*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo407 send*94a8112d-48f8-4239-abdd-5f6ba3d633bd*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo411 send*85575f0a-07bd-47df-848f-2cb3fad91896*0*1*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo415 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
