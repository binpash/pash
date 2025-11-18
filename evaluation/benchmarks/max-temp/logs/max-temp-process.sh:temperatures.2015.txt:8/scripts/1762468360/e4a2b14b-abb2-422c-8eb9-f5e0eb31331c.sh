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
{ rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo541" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo542" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo545" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo546" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo549" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo550" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo553" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo554" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo557" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo558" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo561" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo562" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo659" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo663" ; } 
 { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo667" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo541" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo542" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo545" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo546" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo549" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo550" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo553" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo554" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo557" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo558" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo561" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo562" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo659" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo663" ; } 
 { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo667" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ bigram_aux_reduce "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo542" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo546" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo550" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo554" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo558" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo562" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo659" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo663" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo667" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo541" -o "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo542" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo545" -o "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo546" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo549" -o "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo550" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo553" -o "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo554" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo557" -o "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo558" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo561" -o "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo562" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*ecc3bc6b-a1e2-44ce-ad7c-a9fe90849fcf*1*0*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo541 recv*ca7c9f3b-0b1c-444b-a117-bcfdd65fde84*1*0*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo553 send*09eec0cc-f889-4e5a-a09c-d3aa243ea454*0*1*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo659 & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*82204831-dacc-4c34-a3a8-9763b5424ac7*1*0*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo545 recv*6fa151a2-6520-407f-97f8-02b93bec06e4*1*0*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo557 send*5a1bf929-5470-410b-bc30-e2375890e9ba*0*1*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo663 & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*71419a04-25cc-4675-819c-c276bbab3df6*1*0*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo549 recv*b87b73fd-b081-4057-87d1-593e5452ca38*1*0*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo561 send*fdcf033d-291e-4b67-87cd-f147c52cb121*0*1*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo667 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
