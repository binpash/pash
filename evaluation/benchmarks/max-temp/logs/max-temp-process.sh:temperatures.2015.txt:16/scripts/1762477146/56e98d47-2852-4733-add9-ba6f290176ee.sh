export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_wHLNLXG/ 
mkdir -p /tmp/pash_wHLNLXG/7e971a3ac9aa4b539c901882d7a153a2/ 
mkdir -p /tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/ 
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
{ rm -f "/tmp/pash_wHLNLXG/7e971a3ac9aa4b539c901882d7a153a2/#fifo2" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo57" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo61" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo65" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo69" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo73" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo77" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo81" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo85" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo89" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo93" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo97" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo101" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo105" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo109" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo113" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo117" ; } 
 { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo188" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_wHLNLXG/7e971a3ac9aa4b539c901882d7a153a2/#fifo2" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo57" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo61" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo65" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo69" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo73" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo77" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo81" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo85" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo89" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo93" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo97" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo101" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo105" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo109" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo113" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo117" ; } 
 { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo188" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ cat "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo188" >"/tmp/pash_wHLNLXG/7e971a3ac9aa4b539c901882d7a153a2/#fifo2" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_split "/tmp/pash_wHLNLXG/7e971a3ac9aa4b539c901882d7a153a2/#fifo2" 1000000 "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo57" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo61" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo65" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo69" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo73" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo77" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo81" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo85" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo89" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo93" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo97" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo101" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo105" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo109" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo113" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo117" & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3.9 aws/s3-get-object.py max-temp/inputs/temperatures.2015.txt "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo188" & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib send*b9001e18-e05d-4dfb-8400-fb15abb2d0f2*0*1*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo57 send*0e89146e-4ecb-4f59-a7b1-fb102c2b520d*0*1*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo61 send*ad11c6b0-f61d-4333-8ba2-77e38f1eac2c*0*1*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo65 send*0cd83453-533d-4902-8cec-0cf7e3b8f211*0*1*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo69 send*50b606fd-6bb7-48c7-b0c6-576869219d72*0*1*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo73 send*16e9f772-3815-4216-8ff6-0ed756f2c882*0*1*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo77 send*9ba3beec-abdf-4cdb-8bbb-f901f0487ac5*0*1*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo81 send*0950e8c7-cb40-4f73-9db6-7c89f544ba9d*0*1*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo85 send*af69c101-85f0-4d97-8755-062197b49a39*0*1*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo89 send*757cfae5-fad7-460c-913a-bce54b33f74c*0*1*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo93 send*8dc99df0-5f6d-4086-902b-d3fc0d0a229d*0*1*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo97 send*8340111f-2ae5-4a1a-8937-05de5a246782*0*1*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo101 send*8dba8f86-f648-495a-91b2-324393569bab*0*1*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo105 send*dba8ba93-417b-4265-a27f-1e2bdca34be8*0*1*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo109 send*745c6270-4eae-4a24-9274-28620eacb478*0*1*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo113 send*8078f800-6ae9-45f8-9b47-dc934b6b9cc1*0*1*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo117 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
