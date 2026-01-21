export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export LOCPATH=/var/task/runtime/locale
export LANG=C.UTF-8
export LOCALE_ALL=C.UTF-8
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_jGtSG4s/ 
mkdir -p /tmp/pash_jGtSG4s/8eb678a9aafb49c0b4d6ea8210fe722c/ 
mkdir -p /tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/ 
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
{ rm -f "/tmp/pash_jGtSG4s/8eb678a9aafb49c0b4d6ea8210fe722c/#fifo20" ; }
 { { rm -f "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo35" ; }
 { { rm -f "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo35_temp" ; }
 { { rm -f "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo59" ; }
 { rm -f /tmp/debug_l2_output.tmp ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_jGtSG4s/8eb678a9aafb49c0b4d6ea8210fe722c/#fifo20" ; }
 { { mkfifo "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo35" ; }
 { { mkfifo "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo35_temp" ; }
 { mkfifo "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo59" ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
touch /tmp/debug_l2_output.tmp
pids_to_kill=""
{ bash -c ' sed "s/T..:..:..//" ' <"/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo59" >"/tmp/pash_jGtSG4s/8eb678a9aafb49c0b4d6ea8210fe722c/#fifo20" & }
pids_to_kill="${!} ${pids_to_kill}"
{ bash -c ' cut -d "," -f 1,3 ' <"/tmp/pash_jGtSG4s/8eb678a9aafb49c0b4d6ea8210fe722c/#fifo20" >"/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo35_temp" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo35_temp" -o "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo35" -o /tmp/debug_l2_output.tmp -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3 aws/s3-put-object.py "debug-logs/$version/l2-output.txt" /tmp/debug_l2_output.tmp "$version" & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3 aws/s3-shard-reader.py "covid-mts/inputs/in_tiny.csv" "/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo59" bytes=4000028-8000056 shard=1 num_shards=2 job_uid=j1 debug=True & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib send*bb*0*1*/tmp/pash_jGtSG4s/0da3b8886e4e4202bc76cb60f6f320eb/#fifo35 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )



# optimize
# pattern matching file -> rsplit -> outs
# remove that first subgrpa 
# outs: direct streaming 
# first subgraph - ... - > merger EC2 node : PATH p
# p : remove r_wraps and replace r_merge with cat 
