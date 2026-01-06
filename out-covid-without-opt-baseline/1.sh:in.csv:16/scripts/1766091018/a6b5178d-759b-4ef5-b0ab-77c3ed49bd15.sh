export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_gbcFrw6/ 
mkdir -p /tmp/pash_gbcFrw6/cd4ea6dd729a4293aa848eb873385dd9/ 
mkdir -p /tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/ 
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
{ rm -f "/tmp/pash_gbcFrw6/cd4ea6dd729a4293aa848eb873385dd9/#fifo12" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/cd4ea6dd729a4293aa848eb873385dd9/#fifo14" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo323" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo324" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo327" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo328" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo331" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo332" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo335" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo336" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo339" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo340" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo343" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo344" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo347" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo348" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo351" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo352" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo355" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo356" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo359" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo360" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo363" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo364" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo367" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo368" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo371" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo372" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo375" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo376" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo379" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo380" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo383" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo384" ; } 
 { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo385" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_gbcFrw6/cd4ea6dd729a4293aa848eb873385dd9/#fifo12" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/cd4ea6dd729a4293aa848eb873385dd9/#fifo14" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo323" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo324" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo327" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo328" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo331" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo332" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo335" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo336" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo339" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo340" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo343" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo344" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo347" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo348" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo351" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo352" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo355" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo356" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo359" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo360" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo363" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo364" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo367" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo368" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo371" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo372" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo375" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo376" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo379" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo380" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo383" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo384" ; } 
 { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo385" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ sort -m "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo324" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo328" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo332" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo336" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo340" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo344" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo348" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo352" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo356" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo360" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo364" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo368" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo372" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo376" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo380" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo384" >"/tmp/pash_gbcFrw6/cd4ea6dd729a4293aa848eb873385dd9/#fifo12" & }
pids_to_kill="${!} ${pids_to_kill}"
{ uniq -c <"/tmp/pash_gbcFrw6/cd4ea6dd729a4293aa848eb873385dd9/#fifo12" >"/tmp/pash_gbcFrw6/cd4ea6dd729a4293aa848eb873385dd9/#fifo14" & }
pids_to_kill="${!} ${pids_to_kill}"
{ awk "{print \$2,\$1}" <"/tmp/pash_gbcFrw6/cd4ea6dd729a4293aa848eb873385dd9/#fifo14" >"/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo385" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo323" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo324" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo327" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo328" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo331" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo332" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo335" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo336" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo339" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo340" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo343" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo344" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo347" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo348" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo351" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo352" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo355" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo356" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo359" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo360" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo363" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo364" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo367" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo368" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo371" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo372" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo375" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo376" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo379" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo380" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo383" -o "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo384" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*af79f8f0-3540-4821-866a-f9509d679833*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo323 recv*0a42fa2f-335e-4663-bab0-1e7a8a22c3c7*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo327 recv*4c9b170f-f98a-4979-bd9a-b3b2a297c499*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo331 recv*816bf289-d101-4a3f-b1a2-e04e4943e857*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo335 recv*61d8554c-731a-43cc-a921-aba86a5afdb5*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo339 recv*eb49f9c6-e509-4885-831e-13ec092f66d8*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo343 recv*ba0cf58d-fad2-4fc7-9419-1669a13a5b93*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo347 recv*a441cf33-8f32-47ad-a27e-060a306e7889*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo351 recv*5204950d-b55a-400c-adc4-95902894b64b*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo355 recv*007bf207-316f-4efa-8abf-66336e68d1a5*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo359 recv*cd2f8cf7-6a66-404b-b9b1-47444ce8653d*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo363 recv*94e445cf-7cac-493b-ac6b-115867e62fa8*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo367 recv*401c72c6-332f-4c61-9cc6-67d00270db59*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo371 recv*95bb98e7-a5f4-487e-8c59-53a5c55a81eb*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo375 recv*02c50ce4-0daa-4b6d-9a81-0f38408911f7*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo379 recv*536ed251-cb9e-4358-8f07-1c1c6dad0616*1*0*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo383 & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3.9 aws/s3-put-object.py covid-mts/outputs/1.sh:in.csv:16:hybridstdout.txt "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo385" $1 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
