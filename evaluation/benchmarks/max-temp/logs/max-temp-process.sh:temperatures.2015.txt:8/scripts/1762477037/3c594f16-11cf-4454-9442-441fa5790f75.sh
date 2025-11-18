export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_1hm9dlJ/ 
mkdir -p /tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/ 
mkdir -p /tmp/pash_1hm9dlJ/794a585ac5a44193a676f633c3bf21a7/ 
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
{ rm -f "/tmp/pash_1hm9dlJ/794a585ac5a44193a676f633c3bf21a7/#fifo2" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo58" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo62" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo66" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo70" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo74" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo78" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo82" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo86" ; } 
 { rm -f "/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo149" ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_1hm9dlJ/794a585ac5a44193a676f633c3bf21a7/#fifo2" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo58" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo62" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo66" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo70" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo74" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo78" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo82" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo86" ; } 
 { mkfifo "/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo149" ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ cat "/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo149" >"/tmp/pash_1hm9dlJ/794a585ac5a44193a676f633c3bf21a7/#fifo2" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_split "/tmp/pash_1hm9dlJ/794a585ac5a44193a676f633c3bf21a7/#fifo2" 1000000 "/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo58" "/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo62" "/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo66" "/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo70" "/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo74" "/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo78" "/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo82" "/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo86" & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3.9 aws/s3-get-object.py max-temp/inputs/temperatures.2015.txt "/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo149" & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib send*2e762e31-a71e-4d6c-b350-3f7538e9e2fc*0*1*/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo58 send*6f2282dd-1068-4127-860e-c16e07e7e224*0*1*/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo62 send*8765bf68-362b-464a-bce6-1e4527655ae9*0*1*/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo66 send*8b272cc1-b7a4-4841-91f3-1166561fb472*0*1*/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo70 send*01e47a5a-0916-4a48-9ef5-12d46c2de61e*0*1*/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo74 send*0acc82ef-9f89-4e67-995b-72d9d7249c7b*0*1*/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo78 send*5bc7e00c-4774-4a27-83d5-fe11c20aebab*0*1*/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo82 send*61e2744f-31c6-4d37-917e-b9e9bf5a5622*0*1*/tmp/pash_1hm9dlJ/167cb7ca2e5c46dab193004b9f65c53d/#fifo86 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
