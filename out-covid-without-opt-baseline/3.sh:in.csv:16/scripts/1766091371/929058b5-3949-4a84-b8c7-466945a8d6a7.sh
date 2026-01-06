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
{ rm -f "/tmp/pash_psR0tRP/8f01866ff5f743a8b65995a89637b3f5/#fifo16" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo485" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo486" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo489" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo490" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo493" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo494" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo497" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo498" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo501" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo502" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo505" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo506" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo509" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo510" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo513" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo514" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo517" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo518" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo521" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo522" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo525" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo526" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo529" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo530" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo533" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo534" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo537" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo538" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo541" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo542" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo545" ; } 
 { { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo546" ; } 
 { rm -f "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo547" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_psR0tRP/8f01866ff5f743a8b65995a89637b3f5/#fifo16" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo485" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo486" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo489" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo490" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo493" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo494" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo497" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo498" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo501" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo502" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo505" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo506" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo509" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo510" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo513" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo514" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo517" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo518" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo521" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo522" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo525" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo526" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo529" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo530" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo533" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo534" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo537" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo538" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo541" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo542" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo545" ; } 
 { { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo546" ; } 
 { mkfifo "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo547" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ sort -k 1 -n -m "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo486" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo490" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo494" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo498" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo502" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo506" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo510" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo514" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo518" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo522" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo526" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo530" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo534" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo538" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo542" "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo546" >"/tmp/pash_psR0tRP/8f01866ff5f743a8b65995a89637b3f5/#fifo16" & }
pids_to_kill="${!} ${pids_to_kill}"
{ awk "{print \$2,\$1}" <"/tmp/pash_psR0tRP/8f01866ff5f743a8b65995a89637b3f5/#fifo16" >"/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo547" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo485" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo486" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo489" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo490" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo493" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo494" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo497" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo498" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo501" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo502" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo505" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo506" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo509" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo510" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo513" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo514" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo517" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo518" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo521" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo522" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo525" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo526" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo529" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo530" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo533" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo534" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo537" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo538" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo541" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo542" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo545" -o "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo546" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*c337e462-e180-4e5d-af4f-0ad0f74fccb4*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo485 recv*b6091cd9-7cf6-4562-9d3f-887ad155b669*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo489 recv*3a971bc5-657a-4edf-b613-02e5676497de*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo493 recv*5bd2e1ab-ac5f-45ee-a2c2-b9b15d517fb4*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo497 recv*480d78cd-5b90-4de3-a0c5-41a285d2f6c3*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo501 recv*9426aed3-235f-4743-b026-9950f1036d9a*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo505 recv*dfafa08e-ec61-47c8-a95f-19fd3ef5d07a*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo509 recv*bc887b1d-e509-495e-bd98-0b046770faaa*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo513 recv*ed2a3589-4f2e-4ab8-b9ae-f315a394747e*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo517 recv*7e959eef-5dc7-459d-8e5e-14a84d75a1dc*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo521 recv*f97cf89f-72cd-41a8-8aa5-8e36e4269cdf*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo525 recv*f92f3a9e-fd78-4494-acd6-b3b5b054a590*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo529 recv*97032130-92ab-484e-9db9-cb11c0960686*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo533 recv*6a351920-aa10-457a-904b-f48757c73ef6*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo537 recv*abba247b-4fe1-4f52-af76-bd299c956cd1*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo541 recv*aa7bf13e-8d46-4365-adf7-c31eae807dff*1*0*/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo545 & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3.9 aws/s3-put-object.py "covid-mts/outputs/3.sh:in.csv:16:hybrid"stdout.txt "/tmp/pash_psR0tRP/67e763494cdc4774bef766b1dde061fa/#fifo547" $1 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
