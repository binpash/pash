export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_1hm9dlJ/ 
mkdir -p /tmp/pash_1hm9dlJ/9fbf32595d454f71adf7a269057d2b79/ 
mkdir -p /tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/ 
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
{ rm -f "/tmp/pash_1hm9dlJ/9fbf32595d454f71adf7a269057d2b79/#fifo2" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo33" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo37" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo41" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo45" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo49" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo53" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo57" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo61" ; } 
 { rm -f "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo100" ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_1hm9dlJ/9fbf32595d454f71adf7a269057d2b79/#fifo2" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo33" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo37" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo41" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo45" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo49" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo53" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo57" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo61" ; } 
 { mkfifo "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo100" ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ cat "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo100" >"/tmp/pash_1hm9dlJ/9fbf32595d454f71adf7a269057d2b79/#fifo2" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_split "/tmp/pash_1hm9dlJ/9fbf32595d454f71adf7a269057d2b79/#fifo2" 1000000 "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo33" "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo37" "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo41" "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo45" "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo49" "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo53" "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo57" "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo61" & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3.9 aws/s3-get-object.py max-temp/inputs/temperatures.2015.txt "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo100" & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib send*adaf05e5-93f0-4f42-9d1b-58d0e09261b5*0*1*/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo33 send*fb9fdba3-39e9-445a-9f63-d8a265d41273*0*1*/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo37 send*73d03806-d797-4d77-a1d8-a5eb10a86ca9*0*1*/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo41 send*93667898-14e6-48a6-9814-418fb99ba8c2*0*1*/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo45 send*0a44ce98-3aa7-45d1-8bf4-6acbf3a05e54*0*1*/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo49 send*3816bda2-233e-4bc9-b042-0dbae245257e*0*1*/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo53 send*c7389dbb-4069-49ab-82a7-fad2db163a49*0*1*/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo57 send*d250b5f5-305d-4b8f-9386-9454c465fb54*0*1*/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo61 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
