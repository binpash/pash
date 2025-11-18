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
{ rm -f "/tmp/pash_wHLNLXG/7e971a3ac9aa4b539c901882d7a153a2/#fifo6" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo123" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo124" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo127" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo128" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo131" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo132" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo135" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo136" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo139" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo140" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo143" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo144" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo147" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo148" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo151" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo152" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo155" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo156" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo159" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo160" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo163" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo164" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo167" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo168" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo171" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo172" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo175" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo176" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo179" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo180" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo183" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo184" ; } 
 { rm -f "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo185" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_wHLNLXG/7e971a3ac9aa4b539c901882d7a153a2/#fifo6" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo123" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo124" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo127" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo128" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo131" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo132" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo135" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo136" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo139" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo140" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo143" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo144" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo147" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo148" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo151" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo152" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo155" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo156" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo159" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo160" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo163" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo164" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo167" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo168" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo171" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo172" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo175" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo176" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo179" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo180" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo183" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo184" ; } 
 { mkfifo "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo185" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_merge "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo124" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo128" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo132" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo136" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo140" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo144" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo148" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo152" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo156" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo160" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo164" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo168" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo172" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo176" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo180" "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo184" >"/tmp/pash_wHLNLXG/7e971a3ac9aa4b539c901882d7a153a2/#fifo6" & }
pids_to_kill="${!} ${pids_to_kill}"
{ awk "{ total += \$1; count++ } END { print total/count }" <"/tmp/pash_wHLNLXG/7e971a3ac9aa4b539c901882d7a153a2/#fifo6" >"/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo185" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo123" -o "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo124" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo127" -o "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo128" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo131" -o "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo132" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo135" -o "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo136" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo139" -o "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo140" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo143" -o "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo144" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo147" -o "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo148" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo151" -o "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo152" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo155" -o "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo156" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo159" -o "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo160" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo163" -o "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo164" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo167" -o "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo168" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo171" -o "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo172" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo175" -o "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo176" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo179" -o "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo180" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo183" -o "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo184" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*53a12b76-b552-492f-b3ec-29dd9a35f9db*1*0*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo123 recv*f38b77ad-0b9d-46b0-b417-b74ef6ef1417*1*0*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo127 recv*31f61e43-b02e-4ee5-ab7d-814d12754960*1*0*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo131 recv*0b860a7f-b04b-4583-a7ca-94a092236b4d*1*0*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo135 recv*dda7c0aa-35d6-481f-a133-ee55cce36b06*1*0*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo139 recv*17cce4c1-bbcd-4b24-b4a4-cddc86f6be37*1*0*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo143 recv*371e89ed-02c2-49f9-9ab9-4b9654550e1e*1*0*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo147 recv*0fff2a3c-8d74-4eaf-95fc-a1400eb327a0*1*0*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo151 recv*7fb1c60f-d2eb-4840-9364-b153c8b45f92*1*0*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo155 recv*53716067-c6d7-42d4-876f-ab57effd4ca5*1*0*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo159 recv*cb98da2a-f3a6-4c7e-a9f2-de893f7a3705*1*0*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo163 recv*1a57ff68-361e-4dcd-9252-8ebaffb7c7eb*1*0*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo167 recv*a5f1a74c-2e32-4ba2-a0f2-31b5a987e350*1*0*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo171 recv*9ded17a2-482e-4510-beb7-66cea15d8b9f*1*0*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo175 recv*69ad6707-5a88-44be-b5fc-03afca800e2a*1*0*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo179 recv*ccf69d6c-6efe-4b81-986e-bbc1876a25e3*1*0*/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo183 & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3.9 aws/s3-put-object.py max-temp/outputs/max-temp-process.sh:temperatures.2015.txt:16:hybrid:average.stdout.txt "/tmp/pash_wHLNLXG/8235ba140a954a29b977780953950992/#fifo185" $1 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
