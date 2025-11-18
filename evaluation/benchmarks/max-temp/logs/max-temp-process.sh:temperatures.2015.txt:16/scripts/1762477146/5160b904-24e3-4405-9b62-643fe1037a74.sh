export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_wHLNLXG/ 
mkdir -p /tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/ 
mkdir -p /tmp/pash_wHLNLXG/cc82ba3224fb4321bfa669ada022e3e7/ 
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
{ rm -f "/tmp/pash_wHLNLXG/cc82ba3224fb4321bfa669ada022e3e7/#fifo2" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo106" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo110" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo114" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo118" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo122" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo126" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo130" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo134" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo138" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo142" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo146" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo150" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo154" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo158" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo162" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo166" ; } 
 { rm -f "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo293" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_wHLNLXG/cc82ba3224fb4321bfa669ada022e3e7/#fifo2" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo106" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo110" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo114" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo118" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo122" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo126" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo130" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo134" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo138" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo142" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo146" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo150" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo154" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo158" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo162" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo166" ; } 
 { mkfifo "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo293" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ cat "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo293" >"/tmp/pash_wHLNLXG/cc82ba3224fb4321bfa669ada022e3e7/#fifo2" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_split "/tmp/pash_wHLNLXG/cc82ba3224fb4321bfa669ada022e3e7/#fifo2" 1000000 "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo106" "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo110" "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo114" "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo118" "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo122" "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo126" "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo130" "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo134" "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo138" "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo142" "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo146" "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo150" "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo154" "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo158" "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo162" "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo166" & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3.9 aws/s3-get-object.py max-temp/inputs/temperatures.2015.txt "/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo293" & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib send*b6789ce5-e915-4773-85bb-cb3c22118f2d*0*1*/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo106 send*e345f497-2710-43bc-8671-3840c066289d*0*1*/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo110 send*11b088d4-50ed-43ab-a449-b68491835548*0*1*/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo114 send*00f30339-79a4-4ef9-9d72-145891d881a9*0*1*/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo118 send*5fcc1b42-8cbb-482e-9b75-356d3fb8c324*0*1*/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo122 send*3443d591-b0a2-4de2-9ad0-35dd3d4409bc*0*1*/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo126 send*51b77278-611f-4604-82da-67bd651e0a01*0*1*/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo130 send*21232e19-3f5c-40a5-bcf7-952d082364be*0*1*/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo134 send*f11107b1-9b0f-487c-a238-410c84236441*0*1*/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo138 send*8a4fba82-61e0-4d0b-bf08-f9f145e220a4*0*1*/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo142 send*f1fe5271-f45f-450a-b007-2d21f920e0c8*0*1*/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo146 send*b292a535-5231-489b-847d-41473d9333b1*0*1*/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo150 send*9ba195d5-a618-4a2b-b130-7777a4067e32*0*1*/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo154 send*ba033f2a-1325-4033-a3cc-add4f47f0a0a*0*1*/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo158 send*18bcfc00-cacc-4cb0-92e6-12ad980ae226*0*1*/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo162 send*aefbf667-4d3e-4778-ac33-5989e5a2c545*0*1*/tmp/pash_wHLNLXG/f67e6301e48b49f1bd3ffa14021a9848/#fifo166 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
