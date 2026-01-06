export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_vtZsjBm/ 
mkdir -p /tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/ 
mkdir -p /tmp/pash_vtZsjBm/68a7f4eeec164679bf1bc1313be933aa/ 
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
{ rm -f "/tmp/pash_vtZsjBm/68a7f4eeec164679bf1bc1313be933aa/#fifo12" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/68a7f4eeec164679bf1bc1313be933aa/#fifo14" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo323" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo324" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo327" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo328" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo331" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo332" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo335" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo336" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo339" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo340" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo343" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo344" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo347" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo348" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo351" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo352" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo355" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo356" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo359" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo360" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo363" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo364" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo367" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo368" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo371" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo372" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo375" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo376" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo379" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo380" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo383" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo384" ; } 
 { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo385" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_vtZsjBm/68a7f4eeec164679bf1bc1313be933aa/#fifo12" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/68a7f4eeec164679bf1bc1313be933aa/#fifo14" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo323" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo324" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo327" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo328" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo331" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo332" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo335" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo336" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo339" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo340" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo343" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo344" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo347" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo348" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo351" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo352" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo355" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo356" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo359" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo360" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo363" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo364" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo367" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo368" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo371" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo372" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo375" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo376" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo379" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo380" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo383" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo384" ; } 
 { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo385" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ sort -m "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo324" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo328" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo332" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo336" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo340" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo344" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo348" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo352" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo356" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo360" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo364" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo368" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo372" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo376" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo380" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo384" >"/tmp/pash_vtZsjBm/68a7f4eeec164679bf1bc1313be933aa/#fifo12" & }
pids_to_kill="${!} ${pids_to_kill}"
{ uniq -c <"/tmp/pash_vtZsjBm/68a7f4eeec164679bf1bc1313be933aa/#fifo12" >"/tmp/pash_vtZsjBm/68a7f4eeec164679bf1bc1313be933aa/#fifo14" & }
pids_to_kill="${!} ${pids_to_kill}"
{ awk "{print \$2,\$1}" <"/tmp/pash_vtZsjBm/68a7f4eeec164679bf1bc1313be933aa/#fifo14" >"/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo385" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo323" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo324" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo327" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo328" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo331" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo332" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo335" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo336" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo339" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo340" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo343" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo344" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo347" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo348" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo351" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo352" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo355" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo356" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo359" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo360" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo363" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo364" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo367" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo368" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo371" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo372" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo375" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo376" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo379" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo380" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo383" -o "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo384" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*2dec585e-ce00-436d-83c2-a99cdaebe001*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo323 recv*a9e06696-8c97-4342-9fdb-df1f3a94ba7f*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo327 recv*dd950298-bfed-48a7-8abb-cca797d06cf5*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo331 recv*d037964c-e3a6-4546-b783-0a8ccb0b7fce*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo335 recv*b69a2b1f-1d84-458b-9b0e-4399ca98fc8d*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo339 recv*5d8c2366-be70-410d-9b5f-a3a2264f51e1*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo343 recv*2a7bd908-6063-4f38-821b-5e13471b6ed9*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo347 recv*261b3b22-0a98-4b81-a51d-847b3e71c1ab*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo351 recv*4c019f8e-3473-497c-b8ef-77a6383353b1*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo355 recv*87e3a5a4-b7f1-4fcd-b1c0-5e1e36990490*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo359 recv*87c6a072-64ca-4877-a625-b1095b88429c*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo363 recv*a600c370-b11c-4cc7-9b58-17b9af73a6c8*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo367 recv*6741e953-6a79-4f9b-96ce-51fe82e6ed4c*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo371 recv*e75bfd00-c047-4702-aa09-0b03c2a54778*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo375 recv*fed3c98b-0c12-48cb-88c6-14d8cdd97524*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo379 recv*d42ade91-e6c0-415e-bbf0-db7989833b3d*1*0*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo383 & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3.9 aws/s3-put-object.py "covid-mts/outputs/4.sh:in.csv:16:hybrid"stdout.txt "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo385" $1 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
