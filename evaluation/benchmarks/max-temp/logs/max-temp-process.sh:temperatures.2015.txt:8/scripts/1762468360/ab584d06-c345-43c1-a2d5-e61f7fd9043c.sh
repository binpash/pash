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
{ rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo397" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo398" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo401" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo402" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo405" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo406" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo409" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo410" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo413" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo414" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo417" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo418" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo587" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo591" ; } 
 { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo595" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo397" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo398" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo401" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo402" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo405" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo406" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo409" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo410" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo413" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo414" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo417" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo418" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo587" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo591" ; } 
 { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo595" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ bigram_aux_reduce "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo398" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo402" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo406" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo410" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo414" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo418" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo587" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo591" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo595" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo397" -o "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo398" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo401" -o "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo402" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo405" -o "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo406" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo409" -o "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo410" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo413" -o "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo414" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo417" -o "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo418" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*a68d6048-3022-45e5-9e6c-ebb0510ec30a*1*0*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo397 recv*b1db0761-6242-45b9-be62-33aa4795936b*1*0*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo409 send*25169184-f879-4b5e-aec0-0668d610954d*0*1*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo587 & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*7f2b891c-6a3f-4747-94a2-6babdc578864*1*0*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo401 recv*e3131660-619e-4381-a77f-a2de3a2ac001*1*0*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo413 send*e3d74d5d-2161-4fa5-98ce-e5dff08d511d*0*1*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo591 & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*ae8a8698-cbc7-4ac8-b407-d17764cb1962*1*0*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo405 recv*54121c27-ef62-4d76-82c2-c930f1fe925f*1*0*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo417 send*d6e3e3b5-96f9-4f42-b79b-bce15301985a*0*1*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo595 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
