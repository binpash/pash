export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_cahH684/ 
mkdir -p /tmp/pash_cahH684/566f3e599dfa4a999b957742bec6928b/ 
mkdir -p /tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/ 
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
{ rm -f "/tmp/pash_cahH684/566f3e599dfa4a999b957742bec6928b/#fifo16" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo485" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo486" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo489" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo490" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo493" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo494" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo497" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo498" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo501" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo502" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo505" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo506" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo509" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo510" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo513" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo514" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo517" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo518" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo521" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo522" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo525" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo526" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo529" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo530" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo533" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo534" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo537" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo538" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo541" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo542" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo545" ; } 
 { { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo546" ; } 
 { rm -f "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo547" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_cahH684/566f3e599dfa4a999b957742bec6928b/#fifo16" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo485" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo486" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo489" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo490" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo493" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo494" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo497" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo498" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo501" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo502" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo505" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo506" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo509" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo510" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo513" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo514" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo517" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo518" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo521" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo522" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo525" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo526" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo529" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo530" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo533" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo534" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo537" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo538" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo541" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo542" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo545" ; } 
 { { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo546" ; } 
 { mkfifo "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo547" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ sort -k 1 -n -m "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo486" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo490" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo494" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo498" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo502" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo506" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo510" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo514" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo518" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo522" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo526" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo530" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo534" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo538" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo542" "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo546" >"/tmp/pash_cahH684/566f3e599dfa4a999b957742bec6928b/#fifo16" & }
pids_to_kill="${!} ${pids_to_kill}"
{ awk "{print \$2,\$1}" <"/tmp/pash_cahH684/566f3e599dfa4a999b957742bec6928b/#fifo16" >"/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo547" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo485" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo486" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo489" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo490" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo493" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo494" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo497" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo498" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo501" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo502" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo505" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo506" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo509" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo510" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo513" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo514" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo517" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo518" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo521" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo522" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo525" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo526" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo529" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo530" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo533" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo534" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo537" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo538" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo541" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo542" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo545" -o "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo546" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*4e7dad16-4e28-4949-8c46-0a70d70fd160*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo485 recv*c7a360e8-8bff-4c76-9794-2a2360f1f4f3*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo489 recv*12288f2e-f0a5-4661-a771-72d78f24517a*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo493 recv*aacd2f68-7a0b-40fc-abee-a226153de89f*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo497 recv*707a29e6-973a-4b9e-bb6f-bd5e2ea32f97*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo501 recv*bd56ad31-bfd1-440d-a6e5-101671ec767f*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo505 recv*7b9d1037-12e4-45a3-b765-742cdcdaaa6a*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo509 recv*f75f96b5-e24a-46fd-a8df-2c6e0ec0c707*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo513 recv*d914a39c-318b-4ae1-ba15-4f2d80eafe9c*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo517 recv*a8d36d87-b51e-4a65-9aba-84e399afcde3*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo521 recv*361e908e-27e1-4fc8-a2b1-1112eedffb25*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo525 recv*7603e352-237a-4dbb-88fd-cf007f22fbd6*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo529 recv*b5cbdcba-c2e7-4a83-b5b7-a4cffdc2de65*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo533 recv*77818256-9e8b-4a32-b990-9c84e2441228*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo537 recv*1cf93053-7c07-4fd9-ad15-07d25a87839b*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo541 recv*4e0309a1-10ca-4d42-a29c-cdf1e2c9d86b*1*0*/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo545 & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3.9 aws/s3-put-object.py "covid-mts/outputs/2.sh:in.csv:16:hybrid"stdout.txt "/tmp/pash_cahH684/84c2c7eff0724af883a6c81fb90cbce2/#fifo547" $1 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
