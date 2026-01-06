export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_qrxfus7/ 
mkdir -p /tmp/pash_qrxfus7/9f66d8e02d6a4813aed714a670e9f9a1/ 
mkdir -p /tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/ 
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
{ rm -f "/tmp/pash_qrxfus7/9f66d8e02d6a4813aed714a670e9f9a1/#fifo12" ; } 
 { { rm -f "/tmp/pash_qrxfus7/9f66d8e02d6a4813aed714a670e9f9a1/#fifo14" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo259" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo260" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo263" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo264" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo267" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo268" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo271" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo272" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo275" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo276" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo279" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo280" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo283" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo284" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo287" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo288" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo291" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo292" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo295" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo296" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo299" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo300" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo303" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo304" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo307" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo308" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo311" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo312" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo315" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo316" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo319" ; } 
 { { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo320" ; } 
 { rm -f "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo321" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_qrxfus7/9f66d8e02d6a4813aed714a670e9f9a1/#fifo12" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/9f66d8e02d6a4813aed714a670e9f9a1/#fifo14" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo259" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo260" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo263" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo264" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo267" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo268" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo271" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo272" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo275" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo276" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo279" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo280" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo283" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo284" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo287" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo288" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo291" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo292" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo295" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo296" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo299" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo300" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo303" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo304" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo307" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo308" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo311" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo312" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo315" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo316" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo319" ; } 
 { { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo320" ; } 
 { mkfifo "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo321" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ sort -m "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo260" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo264" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo268" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo272" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo276" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo280" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo284" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo288" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo292" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo296" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo300" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo304" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo308" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo312" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo316" "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo320" >"/tmp/pash_qrxfus7/9f66d8e02d6a4813aed714a670e9f9a1/#fifo12" & }
pids_to_kill="${!} ${pids_to_kill}"
{ uniq -c <"/tmp/pash_qrxfus7/9f66d8e02d6a4813aed714a670e9f9a1/#fifo12" >"/tmp/pash_qrxfus7/9f66d8e02d6a4813aed714a670e9f9a1/#fifo14" & }
pids_to_kill="${!} ${pids_to_kill}"
{ awk "{print \$2,\$1}" <"/tmp/pash_qrxfus7/9f66d8e02d6a4813aed714a670e9f9a1/#fifo14" >"/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo321" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo259" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo260" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo263" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo264" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo267" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo268" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo271" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo272" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo275" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo276" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo279" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo280" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo283" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo284" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo287" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo288" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo291" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo292" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo295" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo296" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo299" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo300" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo303" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo304" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo307" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo308" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo311" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo312" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo315" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo316" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo319" -o "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo320" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*46806224-0496-438c-9543-f81b98c12df3*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo259 recv*f5f4a5f3-25e0-4055-b19f-8a94dea2beef*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo263 recv*59391b99-0028-4798-aa42-b76fd36f9ca5*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo267 recv*509135c7-be9b-4d85-bc26-2637b4236953*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo271 recv*353156a8-c15a-4e41-b5a8-daee9a9976dc*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo275 recv*b7c71bf8-459a-4a89-8970-dc1a54ca4aad*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo279 recv*73994eab-efda-484b-96d7-56331ba39cb8*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo283 recv*6a147068-7f72-4d63-ba39-9d09e610af6f*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo287 recv*943efbf0-6f24-4b9d-853d-d59f208d3448*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo291 recv*60f0cd2d-d673-4e1c-b9ab-2b7326900e1f*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo295 recv*385dd48b-1c4f-4d4e-8d3e-53797810fd52*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo299 recv*28e011e5-a95a-40ec-afdb-9f93e9fd4fa9*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo303 recv*ab1efabd-dd03-48c6-9f68-e59e72a02181*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo307 recv*c016a02f-fe40-4a60-b3d3-9c1c997e6ef5*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo311 recv*545dc060-3950-4291-b8a3-09cd532c6a99*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo315 recv*6ef372a6-e99e-4c73-9375-5dffdea7456a*1*0*/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo319 & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3.9 aws/s3-put-object.py "covid-mts/outputs/4.sh:in.csv:16:hybrid"stdout.txt "/tmp/pash_qrxfus7/d2193d13c92247629a2720bacc324110/#fifo321" $1 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
