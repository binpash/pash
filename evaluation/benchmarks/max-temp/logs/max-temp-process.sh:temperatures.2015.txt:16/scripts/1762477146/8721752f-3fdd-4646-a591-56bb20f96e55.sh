export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_wHLNLXG/ 
mkdir -p /tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/ 
mkdir -p /tmp/pash_wHLNLXG/763e199f4a1f445c913aa3c5ab28edd2/ 
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
{ rm -f "/tmp/pash_wHLNLXG/763e199f4a1f445c913aa3c5ab28edd2/#fifo2" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo106" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo110" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo114" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo118" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo122" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo126" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo130" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo134" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo138" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo142" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo146" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo150" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo154" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo158" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo162" ; } 
 { { rm -f "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo166" ; } 
 { rm -f "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo293" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_wHLNLXG/763e199f4a1f445c913aa3c5ab28edd2/#fifo2" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo106" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo110" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo114" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo118" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo122" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo126" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo130" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo134" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo138" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo142" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo146" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo150" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo154" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo158" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo162" ; } 
 { { mkfifo "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo166" ; } 
 { mkfifo "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo293" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ cat "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo293" >"/tmp/pash_wHLNLXG/763e199f4a1f445c913aa3c5ab28edd2/#fifo2" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_split "/tmp/pash_wHLNLXG/763e199f4a1f445c913aa3c5ab28edd2/#fifo2" 1000000 "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo106" "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo110" "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo114" "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo118" "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo122" "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo126" "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo130" "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo134" "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo138" "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo142" "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo146" "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo150" "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo154" "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo158" "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo162" "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo166" & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3.9 aws/s3-get-object.py max-temp/inputs/temperatures.2015.txt "/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo293" & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib send*35a81340-663a-45f6-a70c-f133a227a571*0*1*/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo106 send*48bfce27-71ce-44cb-a1cd-52fbb0e95057*0*1*/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo110 send*f0f69f4a-4676-4944-97ec-efe8e99825e7*0*1*/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo114 send*afdf6e82-5c01-4328-9b2c-c3b01dfc4c8c*0*1*/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo118 send*d6cfe46d-5f2f-4f9a-98dc-8f137ff84a8a*0*1*/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo122 send*c6b29399-5261-4e94-9976-3d99fa29ef22*0*1*/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo126 send*8c908d2f-b30b-4929-a3d2-673356a0ec0d*0*1*/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo130 send*d72d9108-355b-4afd-9dd3-2b1bbadd55f9*0*1*/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo134 send*18a98229-67a3-438d-873e-bdd296a36eaf*0*1*/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo138 send*fdf7fc88-50a2-4dcb-8bd8-e034edb36f0d*0*1*/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo142 send*f3b7ff59-5431-4381-b16b-94a000a6a79c*0*1*/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo146 send*821b7bc2-5011-4acc-939b-14608fafa4ee*0*1*/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo150 send*c68f797f-b097-465c-967e-a103fc0439fd*0*1*/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo154 send*0bcad73e-ac84-463a-8ee1-be1517b38aab*0*1*/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo158 send*916f95c2-83d4-4c76-bec6-ed51305932d9*0*1*/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo162 send*ca5a5f54-da93-4dc6-85ad-49b0ddfe5d29*0*1*/tmp/pash_wHLNLXG/b6c9fe3245a442f1ac83a63b7d8db664/#fifo166 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
