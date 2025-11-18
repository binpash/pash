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
{ rm -f "/tmp/pash_1hm9dlJ/9fbf32595d454f71adf7a269057d2b79/#fifo6" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo67" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo68" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo71" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo72" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo75" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo76" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo79" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo80" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo83" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo84" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo87" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo88" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo91" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo92" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo95" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo96" ; } 
 { rm -f "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo97" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_1hm9dlJ/9fbf32595d454f71adf7a269057d2b79/#fifo6" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo67" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo68" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo71" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo72" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo75" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo76" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo79" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo80" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo83" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo84" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo87" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo88" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo91" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo92" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo95" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo96" ; } 
 { mkfifo "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo97" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_merge "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo68" "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo72" "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo76" "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo80" "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo84" "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo88" "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo92" "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo96" >"/tmp/pash_1hm9dlJ/9fbf32595d454f71adf7a269057d2b79/#fifo6" & }
pids_to_kill="${!} ${pids_to_kill}"
{ awk "{ total += \$1; count++ } END { print total/count }" <"/tmp/pash_1hm9dlJ/9fbf32595d454f71adf7a269057d2b79/#fifo6" >"/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo97" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo67" -o "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo68" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo71" -o "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo72" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo75" -o "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo76" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo79" -o "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo80" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo83" -o "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo84" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo87" -o "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo88" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo91" -o "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo92" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo95" -o "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo96" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*55332d94-75e0-40a6-a50f-923b1436ddb5*1*0*/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo67 recv*f7d0f20f-8843-4a3b-98f2-6df8e59f67bf*1*0*/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo71 recv*0eb91143-8f0c-41c3-b1cc-81e3f9baa09c*1*0*/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo75 recv*fc4e3f7d-56a0-442a-b219-697b84b8f4fc*1*0*/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo79 recv*3569636e-d669-4a8b-96d6-e7f0e891b7f7*1*0*/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo83 recv*7a7d18ca-5494-4994-904d-85c7c04b0cda*1*0*/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo87 recv*b1e5bf65-dca4-433e-99b9-b5ef30d450bc*1*0*/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo91 recv*032799c9-0784-4189-8b38-8b018463adaf*1*0*/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo95 & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3.9 aws/s3-put-object.py max-temp/outputs/max-temp-process.sh:temperatures.2015.txt:8:hybrid:average.stdout.txt "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo97" $1 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
