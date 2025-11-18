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
{ rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo709" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo710" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo713" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo714" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo717" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo718" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo721" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo722" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo725" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo726" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo729" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo730" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo743" ; } 
 { { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo747" ; } 
 { rm -f "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo751" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo709" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo710" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo713" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo714" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo717" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo718" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo721" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo722" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo725" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo726" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo729" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo730" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo743" ; } 
 { { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo747" ; } 
 { mkfifo "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo751" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ bigram_aux_reduce "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo710" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo714" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo718" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo722" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo726" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo730" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo743" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo747" "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo751" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo709" -o "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo710" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo713" -o "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo714" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo717" -o "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo718" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo721" -o "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo722" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo725" -o "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo726" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo729" -o "/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo730" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*a8f18b24-7307-439c-87e2-c5ab02f95bb0*1*0*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo709 recv*213bce23-afc4-4f31-936b-86cbc91e58f0*1*0*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo721 send*feafc463-ed02-48eb-833c-9d65c9ac58ec*0*1*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo743 & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*6ecb57d8-5669-4265-958f-a519591f3cb1*1*0*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo713 recv*14401aea-5f9d-46bc-a89d-15563375f500*1*0*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo725 send*ceb5f9a9-8515-491b-ab66-ac1fc622ff23*0*1*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo747 & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*bb5510c6-4536-42e0-8842-dc8d0e881b80*1*0*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo717 recv*a8c705aa-96ab-4bd4-9940-dafdf0053f3e*1*0*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo729 send*68c2e8f7-2f70-45be-8aab-af45e550143a*0*1*/tmp/pash_ekD6Q5i/006f8626e71a4b69ba8a06639eaa9fa3/#fifo751 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
