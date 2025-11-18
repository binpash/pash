export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_ekD6Q5i/ 
mkdir -p /tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/ 
bigram_aux_map () 
{ 
    IN=$1;
    OUT=$2;
    AUX_HEAD=$3;
    AUX_TAIL=$4;
    s2=$(mktemp -u);
    aux1=$(mktemp -u);
    aux2=$(mktemp -u);
    aux3=$(mktemp -u);
    temp=$(mktemp -u);
    mkfifo "$s2";
    mkfifo "$aux1";
    mkfifo "$aux2";
    mkfifo "$aux3";
    cat "$IN" > "$temp";
    sed '$d' "$temp" > "$aux3" & cat "$temp" | head -n 1 > "$AUX_HEAD" & cat "$temp" | tail -n 1 > "$AUX_TAIL" & cat "$temp" | tail -n +2 | paste "$aux3" - > "$OUT" & wait;
    rm "$temp";
    rm "$s2";
    rm "$aux1";
    rm "$aux2";
    rm "$aux3"
}
declare -fx bigram_aux_map
bigram_aux_reduce () 
{ 
    IN1=$1;
    AUX_HEAD1=$2;
    AUX_TAIL1=$3;
    IN2=$4;
    AUX_HEAD2=$5;
    AUX_TAIL2=$6;
    OUT=$7;
    AUX_HEAD_OUT=$8;
    AUX_TAIL_OUT=$9;
    temp=$(mktemp -u);
    mkfifo "$temp";
    cat "$AUX_HEAD1" > "$AUX_HEAD_OUT" & cat "$AUX_TAIL2" > "$AUX_TAIL_OUT" & paste "$AUX_TAIL1" "$AUX_HEAD2" > "$temp" & cat "$IN1" "$temp" "$IN2" > "$OUT" & wait;
    rm "$temp"
}
declare -fx bigram_aux_reduce
bigrams_aux () 
{ 
    s2=$(mktemp -u);
    mkfifo "$s2";
    tee "$s2" | tail -n +2 | paste "$s2" - | sed '$d';
    rm "$s2"
}
declare -fx bigrams_aux
inform_daemon_exit () 
{ 
    msg="Exit:${process_id}";
    daemon_response=$(pash_communicate_daemon_just_send "$msg")
}
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
run_parallel () 
{ 
    trap inform_daemon_exit SIGTERM SIGINT EXIT;
    export SCRIPT_TO_EXECUTE="$pash_script_to_execute";
    source "$RUNTIME_DIR/pash_restore_state_and_execute.sh"
}

rm_pash_fifos() {
{ rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo493" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo494" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo497" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo498" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo501" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo502" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo505" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo506" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo509" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo510" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo513" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo514" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo635" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo639" ; } 
 { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo643" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo493" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo494" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo497" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo498" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo501" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo502" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo505" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo506" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo509" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo510" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo513" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo514" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo635" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo639" ; } 
 { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo643" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ bigram_aux_reduce "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo494" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo498" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo502" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo506" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo510" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo514" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo635" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo639" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo643" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo493" -o "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo494" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo497" -o "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo498" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo501" -o "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo502" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo505" -o "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo506" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo509" -o "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo510" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo513" -o "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo514" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*f76dfd14-723e-431b-9eb1-2dc54757b27a*1*0*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo493 recv*c12eb481-eec9-4533-8cec-3b73dad0ac46*1*0*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo505 send*894805a2-6f3c-4579-be81-0d2b342c8341*0*1*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo635 & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*f4e38bf8-5f56-4949-9724-eb1def9aacc3*1*0*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo497 recv*60c333c5-a5f3-4909-9cfe-73615f57b853*1*0*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo509 send*6b8f4dfa-8238-48d9-8fe2-14910b791a4b*0*1*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo639 & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*4911e835-ca78-4c7c-ba44-54fc9bdbe6ac*1*0*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo501 recv*29cfdf17-d594-460d-b1e5-949ceab64f58*1*0*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo513 send*10be3285-ee58-44bb-bf6f-657c7f16a626*0*1*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo643 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
