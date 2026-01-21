export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export LOCPATH=/var/task/runtime/locale
export LANG=C.UTF-8
export LOCALE_ALL=C.UTF-8
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_jGtSG4s/ 
mkdir -p /tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/ 
mkdir -p /tmp/pash_jGtSG4s/e394322ca4d7485a9fb4796fb8f03eb3/ 
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
{ rm -f "/tmp/pash_jGtSG4s/e394322ca4d7485a9fb4796fb8f03eb3/#fifo12" ; }
 { { rm -f "/tmp/pash_jGtSG4s/e394322ca4d7485a9fb4796fb8f03eb3/#fifo14" ; }
 { { rm -f "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo49" ; }
 { { rm -f "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo50" ; }
 { { rm -f "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo53" ; }
 { { rm -f "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo54" ; }
 { { rm -f "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo55" ; }
 { { rm -f /tmp/debug_sortmerge_input_sortl4.tmp ; }
 { { rm -f /tmp/debug_sortmerge_input_sortl3.tmp ; }
 { rm -f /tmp/debug_sortmerge_final.tmp ; } ; } ; } ; } ; } ; } ; } ; }
}
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_jGtSG4s/e394322ca4d7485a9fb4796fb8f03eb3/#fifo12" ; }
 { { mkfifo "/tmp/pash_jGtSG4s/e394322ca4d7485a9fb4796fb8f03eb3/#fifo14" ; }
 { { mkfifo "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo49" ; }
 { { mkfifo "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo50" ; }
 { { mkfifo "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo53" ; }
 { { mkfifo "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo54" ; }
 { mkfifo "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo55" ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
touch /tmp/debug_sortmerge_input_sortl3.tmp /tmp/debug_sortmerge_input_sortl4.tmp /tmp/debug_sortmerge_final.tmp
pids_to_kill=""
{ sort -m "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo50" "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo54" >"/tmp/pash_jGtSG4s/e394322ca4d7485a9fb4796fb8f03eb3/#fifo12" & }
pids_to_kill="${!} ${pids_to_kill}"
{ uniq -c <"/tmp/pash_jGtSG4s/e394322ca4d7485a9fb4796fb8f03eb3/#fifo12" | tee /tmp/debug_sortmerge_final.tmp >"/tmp/pash_jGtSG4s/e394322ca4d7485a9fb4796fb8f03eb3/#fifo14" & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3 aws/s3-put-object.py "debug-logs/$version/sortmerge-final.txt" /tmp/debug_sortmerge_final.tmp "$version" & }
pids_to_kill="${!} ${pids_to_kill}"
{ awk "{print \$2,\$1}" <"/tmp/pash_jGtSG4s/e394322ca4d7485a9fb4796fb8f03eb3/#fifo14" >"/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo55" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo49" -o "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo50" -o /tmp/debug_sortmerge_input_sortl4.tmp -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3 aws/s3-put-object.py "debug-logs/$version/sortmerge-input-from-sortl4.txt" /tmp/debug_sortmerge_input_sortl4.tmp "$version" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo53" -o "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo54" -o /tmp/debug_sortmerge_input_sortl3.tmp -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3 aws/s3-put-object.py "debug-logs/$version/sortmerge-input-from-sortl3.txt" /tmp/debug_sortmerge_input_sortl3.tmp "$version" & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*ff*1*0*/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo49 recv*dd*1*0*/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo53 & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3.9 aws/s3-put-object.py covid-mts/outputs/1.sh:in_tiny.csv:2:hybridstdout.txt "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo55" $1 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
