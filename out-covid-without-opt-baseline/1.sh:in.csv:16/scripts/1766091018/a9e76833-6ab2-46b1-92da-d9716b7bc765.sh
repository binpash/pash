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
{ rm -f "/tmp/pash_gbcFrw6/cd4ea6dd729a4293aa848eb873385dd9/#fifo2" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo129" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo133" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo137" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo141" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo145" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo149" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo153" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo157" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo161" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo165" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo169" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo173" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo177" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo181" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo185" ; } 
 { { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo189" ; } 
 { rm -f "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo388" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_gbcFrw6/cd4ea6dd729a4293aa848eb873385dd9/#fifo2" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo129" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo133" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo137" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo141" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo145" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo149" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo153" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo157" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo161" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo165" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo169" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo173" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo177" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo181" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo185" ; } 
 { { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo189" ; } 
 { mkfifo "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo388" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ cat "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo388" >"/tmp/pash_gbcFrw6/cd4ea6dd729a4293aa848eb873385dd9/#fifo2" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_split "/tmp/pash_gbcFrw6/cd4ea6dd729a4293aa848eb873385dd9/#fifo2" 1000000 "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo129" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo133" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo137" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo141" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo145" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo149" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo153" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo157" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo161" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo165" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo169" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo173" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo177" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo181" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo185" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo189" & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3.9 aws/s3-get-object.py "covid-mts/inputs/in.csv" "/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo388" & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib send*b7ea9b49-ba15-4233-ac82-2375f05cc921*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo129 send*4e99aae8-e06f-4557-a1ee-cc724fdf402d*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo133 send*4e6ec424-e827-4151-958d-94b338872bd1*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo137 send*4eae6ba6-89e7-43cf-8b98-19ff32f67a17*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo141 send*9ee4db70-d246-4306-841f-2f5c61285507*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo145 send*0f991212-0182-4014-83f8-274aca360e82*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo149 send*78716180-3cc2-4371-be15-df1de8a01786*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo153 send*7c789264-dd0d-4f67-934a-6f5a9269536d*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo157 send*28f0fa4f-9544-4309-8340-dedadd2857ae*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo161 send*741b91f3-2987-4c5d-934b-6200a2a48f1d*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo165 send*ecbb1d1a-c2c2-4cf2-a98c-8a4f9a48c0d7*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo169 send*2a79434d-c618-4a1d-9a94-5d50013ef3b5*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo173 send*1352d15e-b19d-4200-b950-66efb1802a69*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo177 send*28cc1916-69bc-41d4-bcd1-6c999dd0dac4*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo181 send*fb8290a4-c571-42b8-ab9e-5d5f20a3ad0c*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo185 send*7d3b911a-8192-4f49-b3e2-5671a12267bf*0*1*/tmp/pash_gbcFrw6/ff29e6ab15d8428d86ffd000a9e6ee48/#fifo189 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
