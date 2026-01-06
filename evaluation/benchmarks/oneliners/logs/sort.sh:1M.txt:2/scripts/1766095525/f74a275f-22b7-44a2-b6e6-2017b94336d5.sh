export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_Otz3ZMa/ 
mkdir -p /tmp/pash_Otz3ZMa/a5d65e2258094baf847baabb3eb6002b/ 
mkdir -p /tmp/pash_Otz3ZMa/efb42938ad50474f9c9b1f81c5b41347/ 
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
{ rm -f "/tmp/pash_Otz3ZMa/a5d65e2258094baf847baabb3eb6002b/#fifo20" ; } 
 { { rm -f "/tmp/pash_Otz3ZMa/efb42938ad50474f9c9b1f81c5b41347/#fifo35" ; } 
 { rm -f "/tmp/pash_Otz3ZMa/efb42938ad50474f9c9b1f81c5b41347/#fifo59" ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_Otz3ZMa/a5d65e2258094baf847baabb3eb6002b/#fifo20" ; } 
 { { mkfifo "/tmp/pash_Otz3ZMa/efb42938ad50474f9c9b1f81c5b41347/#fifo35" ; } 
 { mkfifo "/tmp/pash_Otz3ZMa/efb42938ad50474f9c9b1f81c5b41347/#fifo59" ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_wrap bash -c ' sed "s/T..:..:..//" ' <"/tmp/pash_Otz3ZMa/efb42938ad50474f9c9b1f81c5b41347/#fifo59" >"/tmp/pash_Otz3ZMa/a5d65e2258094baf847baabb3eb6002b/#fifo20" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_wrap bash -c ' cut -d "," -f 1,3 ' <"/tmp/pash_Otz3ZMa/a5d65e2258094baf847baabb3eb6002b/#fifo20" >"/tmp/pash_Otz3ZMa/efb42938ad50474f9c9b1f81c5b41347/#fifo35" & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3.9 aws/s3-shard-reader.py "covid-mts/inputs/in_tiny.csv" "/tmp/pash_Otz3ZMa/efb42938ad50474f9c9b1f81c5b41347/#fifo59" bytes=4000028-8000056 shard=1 num_shards=2 job_uid=60e8b288-1093-428a-9471-35bea33a6bc5 debug=True & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib send*6c9245f9-e24d-4760-8d18-73ac99fa77c4*0*1*/tmp/pash_Otz3ZMa/efb42938ad50474f9c9b1f81c5b41347/#fifo35 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
