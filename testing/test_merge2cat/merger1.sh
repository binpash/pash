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
{ rm -f "/tmp/pash_jGtSG4s/e394322ca4d7485a9fb4796fb8f03eb3/#fifo6" ; }
 { { rm -f "/tmp/pash_jGtSG4s/e394322ca4d7485a9fb4796fb8f03eb3/#fifo8" ; }
 { { rm -f "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo33" ; }
 { { rm -f "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo34" ; }
 { { rm -f "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo37" ; }
 { { rm -f "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo38" ; }
 { { rm -f "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo39" ; }
 { { rm -f "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo43" ; }
 { { rm -f /tmp/debug_merger1_input_l1.tmp ; }
 { { rm -f /tmp/debug_merger1_input_l2.tmp ; }
 { rm -f /tmp/debug_merger1_merged.tmp ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_jGtSG4s/e394322ca4d7485a9fb4796fb8f03eb3/#fifo6" ; }
 { { mkfifo "/tmp/pash_jGtSG4s/e394322ca4d7485a9fb4796fb8f03eb3/#fifo8" ; }
 { { mkfifo "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo33" ; }
 { { mkfifo "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo34" ; }
 { { mkfifo "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo37" ; }
 { { mkfifo "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo38" ; }
 { { mkfifo "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo39" ; }
 { mkfifo "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo43" ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
touch /tmp/debug_merger1_input_l1.tmp /tmp/debug_merger1_input_l2.tmp /tmp/debug_merger1_merged.tmp
pids_to_kill=""
{ cat "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo34" "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo38" >"/tmp/pash_jGtSG4s/e394322ca4d7485a9fb4796fb8f03eb3/#fifo6" & }
pids_to_kill="${!} ${pids_to_kill}"
{ sort -u <"/tmp/pash_jGtSG4s/e394322ca4d7485a9fb4796fb8f03eb3/#fifo6" | tee /tmp/debug_merger1_merged.tmp >"/tmp/pash_jGtSG4s/e394322ca4d7485a9fb4796fb8f03eb3/#fifo8" & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3 aws/s3-put-object.py "debug-logs/$version/merger1-merged.txt" /tmp/debug_merger1_merged.tmp "$version" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_split "/tmp/pash_jGtSG4s/e394322ca4d7485a9fb4796fb8f03eb3/#fifo8" 1000000 "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo39" "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo43" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo33" -o "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo34" -o /tmp/debug_merger1_input_l1.tmp -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3 aws/s3-put-object.py "debug-logs/$version/merger1-input-from-l1.txt" /tmp/debug_merger1_input_l1.tmp "$version" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo37" -o "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo38" -o /tmp/debug_merger1_input_l2.tmp -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3 aws/s3-put-object.py "debug-logs/$version/merger1-input-from-l2.txt" /tmp/debug_merger1_input_l2.tmp "$version" & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*aa*1*0*/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo33 recv*bb*1*0*/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo37 send*ee*0*1*/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo39 send*cc*0*1*/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo43 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
