export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_vweH2q0/ 
mkdir -p /tmp/pash_vweH2q0/707a46a3133345c4ab20761398c817f5/ 
mkdir -p /tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/ 
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
{ rm -f "/tmp/pash_vweH2q0/707a46a3133345c4ab20761398c817f5/#fifo12" ; } 
 { { rm -f "/tmp/pash_vweH2q0/707a46a3133345c4ab20761398c817f5/#fifo14" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo259" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo260" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo263" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo264" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo267" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo268" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo271" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo272" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo275" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo276" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo279" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo280" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo283" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo284" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo287" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo288" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo291" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo292" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo295" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo296" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo299" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo300" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo303" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo304" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo307" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo308" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo311" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo312" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo315" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo316" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo319" ; } 
 { { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo320" ; } 
 { rm -f "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo321" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_vweH2q0/707a46a3133345c4ab20761398c817f5/#fifo12" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/707a46a3133345c4ab20761398c817f5/#fifo14" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo259" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo260" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo263" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo264" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo267" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo268" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo271" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo272" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo275" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo276" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo279" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo280" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo283" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo284" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo287" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo288" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo291" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo292" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo295" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo296" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo299" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo300" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo303" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo304" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo307" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo308" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo311" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo312" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo315" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo316" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo319" ; } 
 { { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo320" ; } 
 { mkfifo "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo321" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ sort -m "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo260" "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo264" "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo268" "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo272" "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo276" "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo280" "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo284" "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo288" "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo292" "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo296" "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo300" "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo304" "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo308" "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo312" "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo316" "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo320" >"/tmp/pash_vweH2q0/707a46a3133345c4ab20761398c817f5/#fifo12" & }
pids_to_kill="${!} ${pids_to_kill}"
{ uniq -c <"/tmp/pash_vweH2q0/707a46a3133345c4ab20761398c817f5/#fifo12" >"/tmp/pash_vweH2q0/707a46a3133345c4ab20761398c817f5/#fifo14" & }
pids_to_kill="${!} ${pids_to_kill}"
{ awk "{print \$2,\$1}" <"/tmp/pash_vweH2q0/707a46a3133345c4ab20761398c817f5/#fifo14" >"/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo321" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo259" -o "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo260" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo263" -o "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo264" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo267" -o "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo268" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo271" -o "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo272" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo275" -o "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo276" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo279" -o "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo280" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo283" -o "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo284" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo287" -o "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo288" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo291" -o "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo292" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo295" -o "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo296" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo299" -o "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo300" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo303" -o "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo304" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo307" -o "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo308" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo311" -o "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo312" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo315" -o "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo316" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo319" -o "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo320" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*874b6917-9bfa-4937-96a2-91112cd6db0a*1*0*/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo259 recv*1f3718c1-1c8c-4ecd-9fde-7fd3a892f7a2*1*0*/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo263 recv*ddebb856-62d3-487c-99c5-40bb6eca9d43*1*0*/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo267 recv*a1fbc44b-b143-47ba-8da6-e6b861d250bf*1*0*/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo271 recv*6c2f2db2-e637-480a-a817-a54c6ab1175d*1*0*/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo275 recv*0d1bc755-a8bf-495a-b9c2-4e3b98649f39*1*0*/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo279 recv*c77572a7-07c0-430c-827a-84cbe8a083fb*1*0*/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo283 recv*06b9a351-133f-4df3-9071-4095f36b5789*1*0*/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo287 recv*c8f85023-21d8-4f20-85a3-4f2b7ff477fa*1*0*/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo291 recv*ad26de48-75af-474c-b27e-e205480c8f41*1*0*/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo295 recv*2499ba59-5bcf-4752-90ad-b3f4e75bc82f*1*0*/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo299 recv*9b6a53b1-6ece-47da-948a-614030859fd8*1*0*/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo303 recv*7d89ddca-9198-4caf-b97d-a7545ecf33e0*1*0*/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo307 recv*ee0cd933-260f-478f-ac64-a0d3d865c8b6*1*0*/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo311 recv*05b3cd89-ebe5-4cc5-97fd-55a75010fe3f*1*0*/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo315 recv*ac79a8e3-2a86-4765-bfa2-c04df23b4acd*1*0*/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo319 & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3.9 aws/s3-put-object.py covid-mts/outputs/1.sh:in.csv:16:hybridstdout.txt "/tmp/pash_vweH2q0/1e4dc43c004f4351908a5e1509cdbdea/#fifo321" $1 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
