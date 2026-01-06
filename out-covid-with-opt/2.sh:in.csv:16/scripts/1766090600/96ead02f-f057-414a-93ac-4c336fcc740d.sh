export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_k1vQSpy/ 
mkdir -p /tmp/pash_k1vQSpy/32e98ac9d8d643d6ac226bf864b2a19e/ 
mkdir -p /tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/ 
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
{ rm -f "/tmp/pash_k1vQSpy/32e98ac9d8d643d6ac226bf864b2a19e/#fifo16" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo421" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo422" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo425" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo426" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo429" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo430" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo433" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo434" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo437" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo438" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo441" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo442" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo445" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo446" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo449" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo450" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo453" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo454" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo457" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo458" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo461" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo462" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo465" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo466" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo469" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo470" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo473" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo474" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo477" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo478" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo481" ; } 
 { { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo482" ; } 
 { rm -f "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo483" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_k1vQSpy/32e98ac9d8d643d6ac226bf864b2a19e/#fifo16" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo421" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo422" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo425" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo426" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo429" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo430" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo433" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo434" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo437" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo438" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo441" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo442" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo445" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo446" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo449" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo450" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo453" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo454" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo457" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo458" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo461" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo462" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo465" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo466" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo469" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo470" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo473" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo474" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo477" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo478" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo481" ; } 
 { { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo482" ; } 
 { mkfifo "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo483" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ sort -k 1 -n -m "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo422" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo426" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo430" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo434" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo438" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo442" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo446" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo450" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo454" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo458" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo462" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo466" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo470" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo474" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo478" "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo482" >"/tmp/pash_k1vQSpy/32e98ac9d8d643d6ac226bf864b2a19e/#fifo16" & }
pids_to_kill="${!} ${pids_to_kill}"
{ awk "{print \$2,\$1}" <"/tmp/pash_k1vQSpy/32e98ac9d8d643d6ac226bf864b2a19e/#fifo16" >"/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo483" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo421" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo422" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo425" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo426" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo429" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo430" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo433" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo434" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo437" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo438" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo441" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo442" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo445" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo446" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo449" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo450" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo453" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo454" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo457" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo458" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo461" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo462" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo465" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo466" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo469" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo470" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo473" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo474" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo477" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo478" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo481" -o "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo482" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*d64ae078-98d7-44d9-8440-656a32994159*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo421 recv*29a0742b-3734-414f-8ceb-11e0d468d2dd*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo425 recv*39538d89-e6c3-4257-a2fa-a6aae9c9f46d*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo429 recv*46587f78-5d63-4d44-aabd-b0be284eceda*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo433 recv*408ae60e-b9d8-4ea4-a0ea-311746f716f4*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo437 recv*5a8f428f-c0ed-4d27-bb15-4e3e3e0d2033*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo441 recv*2989bbfb-0f3e-41f9-b221-6291ce7be619*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo445 recv*f0ad86d9-3aba-4182-8d0b-e93df598aa62*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo449 recv*a7a7e2f3-75bf-4826-a4da-80dce493bd21*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo453 recv*32cd5952-9cc2-4283-a9e6-78fe600a4623*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo457 recv*85f42308-cf9c-4da7-ad8e-381a2155a71c*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo461 recv*7362cce4-2b9c-4b0e-9405-ba8a104459c8*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo465 recv*71cf0217-e95b-4a05-928f-503b36d51da4*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo469 recv*69375662-8b3c-4865-91a2-b6738085d01b*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo473 recv*86c5cef7-9b73-449c-a1fc-8b2840cd4194*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo477 recv*6d7eb7df-deb5-4973-acb1-b20ec4d30f06*1*0*/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo481 & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3.9 aws/s3-put-object.py "covid-mts/outputs/2.sh:in.csv:16:hybrid"stdout.txt "/tmp/pash_k1vQSpy/35f7c74aad414bbaafb74ed87c00e1b3/#fifo483" $1 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
