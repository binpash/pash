export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_1hm9dlJ/ 
mkdir -p /tmp/pash_1hm9dlJ/a149ed255fac47448c135ec6eeb9ff57/ 
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
{ rm -f "/tmp/pash_1hm9dlJ/a149ed255fac47448c135ec6eeb9ff57/#fifo20" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo47" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo48" ; } 
 { rm -f "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo77" ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_1hm9dlJ/a149ed255fac47448c135ec6eeb9ff57/#fifo20" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo47" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo48" ; } 
 { mkfifo "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo77" ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_wrap bash -c ' cut -c 88-92 ' <"/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo48" >"/tmp/pash_1hm9dlJ/a149ed255fac47448c135ec6eeb9ff57/#fifo20" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_wrap bash -c ' grep -v 999 ' <"/tmp/pash_1hm9dlJ/a149ed255fac47448c135ec6eeb9ff57/#fifo20" >"/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo77" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo47" -o "/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo48" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*93667898-14e6-48a6-9814-418fb99ba8c2*1*0*/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo47 send*fc4e3f7d-56a0-442a-b219-697b84b8f4fc*0*1*/tmp/pash_1hm9dlJ/a0af71450a234f488bc9aaa4311a7a7d/#fifo77 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
