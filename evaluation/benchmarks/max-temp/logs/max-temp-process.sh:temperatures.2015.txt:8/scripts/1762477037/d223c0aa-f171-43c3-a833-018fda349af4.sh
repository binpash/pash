export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_1hm9dlJ/ 
mkdir -p /tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/ 
mkdir -p /tmp/pash_1hm9dlJ/275b4d0b7e044557b4ed9f5ab2c39b84/ 
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
{ rm -f "/tmp/pash_1hm9dlJ/275b4d0b7e044557b4ed9f5ab2c39b84/#fifo2" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo58" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo62" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo66" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo70" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo74" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo78" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo82" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo86" ; } 
 { rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo149" ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_1hm9dlJ/275b4d0b7e044557b4ed9f5ab2c39b84/#fifo2" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo58" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo62" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo66" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo70" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo74" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo78" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo82" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo86" ; } 
 { mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo149" ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ cat "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo149" >"/tmp/pash_1hm9dlJ/275b4d0b7e044557b4ed9f5ab2c39b84/#fifo2" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_split "/tmp/pash_1hm9dlJ/275b4d0b7e044557b4ed9f5ab2c39b84/#fifo2" 1000000 "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo58" "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo62" "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo66" "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo70" "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo74" "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo78" "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo82" "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo86" & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3.9 aws/s3-get-object.py max-temp/inputs/temperatures.2015.txt "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo149" & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib send*7d737a17-4e6a-4b8c-b698-3fc7caf379f4*0*1*/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo58 send*74fc0e59-6bee-4f4b-a056-2c8698aa9a15*0*1*/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo62 send*5b4f38e0-0459-4e5f-9e51-d7976b8c7ed3*0*1*/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo66 send*6d7b68bd-6710-4413-9787-30c19de963b1*0*1*/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo70 send*1c033c45-b7ca-4337-9c8d-d872a14df338*0*1*/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo74 send*ae43d4bc-728c-4612-b229-6c58cc7dea4f*0*1*/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo78 send*1fb8aef0-db36-4ee6-b082-731dd588510c*0*1*/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo82 send*f5a8fd98-2c08-406f-b273-b92f05cd413d*0*1*/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo86 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
