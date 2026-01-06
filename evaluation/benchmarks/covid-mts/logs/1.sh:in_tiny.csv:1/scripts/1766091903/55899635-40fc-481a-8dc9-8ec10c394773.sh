export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_sWYgFMy/ 
mkdir -p /tmp/pash_sWYgFMy/26d3668d300e4af3bb28878c673d953e/ 
mkdir -p /tmp/pash_sWYgFMy/0b54f1567ba443aa838a2da098862334/ 
mkdir -p /tmp/pash_sWYgFMy/a3c76bafb94d4c5abb9edeace937eaf3/ 
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
{ rm -f "/tmp/pash_sWYgFMy/0b54f1567ba443aa838a2da098862334/#fifo2" ; } 
 { { rm -f "/tmp/pash_sWYgFMy/a3c76bafb94d4c5abb9edeace937eaf3/#fifo17" ; } 
 { { rm -f "/tmp/pash_sWYgFMy/a3c76bafb94d4c5abb9edeace937eaf3/#fifo18" ; } 
 { { rm -f "/tmp/pash_sWYgFMy/a3c76bafb94d4c5abb9edeace937eaf3/#fifo19" ; } 
 { { rm -f "/tmp/pash_sWYgFMy/0b54f1567ba443aa838a2da098862334/#fifo6" ; } 
 { { rm -f "/tmp/pash_sWYgFMy/0b54f1567ba443aa838a2da098862334/#fifo8" ; } 
 { { rm -f "/tmp/pash_sWYgFMy/a3c76bafb94d4c5abb9edeace937eaf3/#fifo20" ; } 
 { { rm -f "/tmp/pash_sWYgFMy/a3c76bafb94d4c5abb9edeace937eaf3/#fifo21" ; } 
 { { rm -f "/tmp/pash_sWYgFMy/a3c76bafb94d4c5abb9edeace937eaf3/#fifo22" ; } 
 { { rm -f "/tmp/pash_sWYgFMy/a3c76bafb94d4c5abb9edeace937eaf3/#fifo23" ; } 
 { { rm -f "/tmp/pash_sWYgFMy/0b54f1567ba443aa838a2da098862334/#fifo12" ; } 
 { { rm -f "/tmp/pash_sWYgFMy/0b54f1567ba443aa838a2da098862334/#fifo14" ; } 
 { { rm -f "/tmp/pash_sWYgFMy/26d3668d300e4af3bb28878c673d953e/#fifo24" ; } 
 { rm -f "/tmp/pash_sWYgFMy/26d3668d300e4af3bb28878c673d953e/#fifo27" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_sWYgFMy/0b54f1567ba443aa838a2da098862334/#fifo2" ; } 
 { { mkfifo "/tmp/pash_sWYgFMy/a3c76bafb94d4c5abb9edeace937eaf3/#fifo17" ; } 
 { { mkfifo "/tmp/pash_sWYgFMy/a3c76bafb94d4c5abb9edeace937eaf3/#fifo18" ; } 
 { { mkfifo "/tmp/pash_sWYgFMy/a3c76bafb94d4c5abb9edeace937eaf3/#fifo19" ; } 
 { { mkfifo "/tmp/pash_sWYgFMy/0b54f1567ba443aa838a2da098862334/#fifo6" ; } 
 { { mkfifo "/tmp/pash_sWYgFMy/0b54f1567ba443aa838a2da098862334/#fifo8" ; } 
 { { mkfifo "/tmp/pash_sWYgFMy/a3c76bafb94d4c5abb9edeace937eaf3/#fifo20" ; } 
 { { mkfifo "/tmp/pash_sWYgFMy/a3c76bafb94d4c5abb9edeace937eaf3/#fifo21" ; } 
 { { mkfifo "/tmp/pash_sWYgFMy/a3c76bafb94d4c5abb9edeace937eaf3/#fifo22" ; } 
 { { mkfifo "/tmp/pash_sWYgFMy/a3c76bafb94d4c5abb9edeace937eaf3/#fifo23" ; } 
 { { mkfifo "/tmp/pash_sWYgFMy/0b54f1567ba443aa838a2da098862334/#fifo12" ; } 
 { { mkfifo "/tmp/pash_sWYgFMy/0b54f1567ba443aa838a2da098862334/#fifo14" ; } 
 { { mkfifo "/tmp/pash_sWYgFMy/26d3668d300e4af3bb28878c673d953e/#fifo24" ; } 
 { mkfifo "/tmp/pash_sWYgFMy/26d3668d300e4af3bb28878c673d953e/#fifo27" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ cat "/tmp/pash_sWYgFMy/26d3668d300e4af3bb28878c673d953e/#fifo27" >"/tmp/pash_sWYgFMy/0b54f1567ba443aa838a2da098862334/#fifo2" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_split "/tmp/pash_sWYgFMy/0b54f1567ba443aa838a2da098862334/#fifo2" 1000000 "/tmp/pash_sWYgFMy/a3c76bafb94d4c5abb9edeace937eaf3/#fifo17" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_wrap bash -c ' sed "s/T..:..:..//" ' <"/tmp/pash_sWYgFMy/a3c76bafb94d4c5abb9edeace937eaf3/#fifo17" >"/tmp/pash_sWYgFMy/a3c76bafb94d4c5abb9edeace937eaf3/#fifo18" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_wrap bash -c ' cut -d "," -f 1,3 ' <"/tmp/pash_sWYgFMy/a3c76bafb94d4c5abb9edeace937eaf3/#fifo18" >"/tmp/pash_sWYgFMy/a3c76bafb94d4c5abb9edeace937eaf3/#fifo19" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_merge "/tmp/pash_sWYgFMy/a3c76bafb94d4c5abb9edeace937eaf3/#fifo19" >"/tmp/pash_sWYgFMy/0b54f1567ba443aa838a2da098862334/#fifo6" & }
pids_to_kill="${!} ${pids_to_kill}"
{ sort -u <"/tmp/pash_sWYgFMy/0b54f1567ba443aa838a2da098862334/#fifo6" >"/tmp/pash_sWYgFMy/0b54f1567ba443aa838a2da098862334/#fifo8" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_split "/tmp/pash_sWYgFMy/0b54f1567ba443aa838a2da098862334/#fifo8" 1000000 "/tmp/pash_sWYgFMy/a3c76bafb94d4c5abb9edeace937eaf3/#fifo20" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_wrap bash -c ' cut -d "," -f 1 ' <"/tmp/pash_sWYgFMy/a3c76bafb94d4c5abb9edeace937eaf3/#fifo20" >"/tmp/pash_sWYgFMy/a3c76bafb94d4c5abb9edeace937eaf3/#fifo21" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_unwrap <"/tmp/pash_sWYgFMy/a3c76bafb94d4c5abb9edeace937eaf3/#fifo21" >"/tmp/pash_sWYgFMy/a3c76bafb94d4c5abb9edeace937eaf3/#fifo22" & }
pids_to_kill="${!} ${pids_to_kill}"
{ sort <"/tmp/pash_sWYgFMy/a3c76bafb94d4c5abb9edeace937eaf3/#fifo22" >"/tmp/pash_sWYgFMy/a3c76bafb94d4c5abb9edeace937eaf3/#fifo23" & }
pids_to_kill="${!} ${pids_to_kill}"
{ sort -m "/tmp/pash_sWYgFMy/a3c76bafb94d4c5abb9edeace937eaf3/#fifo23" >"/tmp/pash_sWYgFMy/0b54f1567ba443aa838a2da098862334/#fifo12" & }
pids_to_kill="${!} ${pids_to_kill}"
{ uniq -c <"/tmp/pash_sWYgFMy/0b54f1567ba443aa838a2da098862334/#fifo12" >"/tmp/pash_sWYgFMy/0b54f1567ba443aa838a2da098862334/#fifo14" & }
pids_to_kill="${!} ${pids_to_kill}"
{ awk "{print \$2,\$1}" <"/tmp/pash_sWYgFMy/0b54f1567ba443aa838a2da098862334/#fifo14" >"/tmp/pash_sWYgFMy/26d3668d300e4af3bb28878c673d953e/#fifo24" & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3.9 aws/s3-get-object.py "covid-mts/inputs/in_tiny.csv" "/tmp/pash_sWYgFMy/26d3668d300e4af3bb28878c673d953e/#fifo27" & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3.9 aws/s3-put-object.py covid-mts/outputs/1.sh:in_tiny.csv:1:hybridstdout.txt "/tmp/pash_sWYgFMy/26d3668d300e4af3bb28878c673d953e/#fifo24" $1 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
