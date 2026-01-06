export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_psR0tRP/ 
mkdir -p /tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/ 
mkdir -p /tmp/pash_psR0tRP/8f01866ff5f743a8b65995a89637b3f5/ 
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

rm_pash_fifos() {
{ rm -f "/tmp/pash_psR0tRP/8f01866ff5f743a8b65995a89637b3f5/#fifo2" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo163" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo167" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo171" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo175" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo179" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo183" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo187" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo191" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo195" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo199" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo203" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo207" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo211" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo215" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo219" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo223" ; } 
 { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo550" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_psR0tRP/8f01866ff5f743a8b65995a89637b3f5/#fifo2" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo163" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo167" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo171" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo175" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo179" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo183" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo187" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo191" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo195" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo199" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo203" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo207" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo211" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo215" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo219" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo223" ; } 
 { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo550" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ cat "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo550" >"/tmp/pash_psR0tRP/8f01866ff5f743a8b65995a89637b3f5/#fifo2" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_split "/tmp/pash_psR0tRP/8f01866ff5f743a8b65995a89637b3f5/#fifo2" 1000000 "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo163" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo167" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo171" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo175" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo179" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo183" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo187" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo191" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo195" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo199" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo203" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo207" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo211" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo215" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo219" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo223" & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3.9 aws/s3-get-object.py "covid-mts/inputs/in.csv" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo550" & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib send*f9833c5f-4190-4957-b0cb-6f09d6be2296*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo163 send*fdc7db17-df0a-4dc9-b117-632f783ec877*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo167 send*63d44230-32ed-47d6-8c0f-22bbdf9147a3*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo171 send*f507c92a-ec23-4c50-b606-f5fee067beb4*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo175 send*6abd09bb-ff99-4071-a2bd-8316705611e7*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo179 send*24152b66-0a96-45fd-809b-3228e2d132f3*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo183 send*a7a2549b-fe8a-41bb-a0e6-f74711873e06*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo187 send*37a40e3c-25b7-448d-9ec0-fb5b52968b8c*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo191 send*2f7f1624-4e8d-4405-9126-857c7948e8b2*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo195 send*7b7ea795-1aea-4460-971a-a26f2cfb31ac*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo199 send*c0627629-d40a-4837-bbde-993d6cc437a3*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo203 send*dd97395e-36f2-4168-9ec3-8c00fee5a497*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo207 send*5aa4f7f4-0cef-4e80-9225-14677ac7fefd*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo211 send*834e07cd-07dc-426e-a83c-e2d03e65bfc1*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo215 send*508fd1d3-dcdf-4eb6-8ace-589d45798aba*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo219 send*2b24a632-e823-4057-aed4-6c7b48eaf6f1*0*1*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo223 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
