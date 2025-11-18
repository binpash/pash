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
{ rm -f "/tmp/pash_1hm9dlJ/275b4d0b7e044557b4ed9f5ab2c39b84/#fifo8" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo140" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo141" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo144" ; } 
 { { rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo145" ; } 
 { rm -f "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo146" ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_1hm9dlJ/275b4d0b7e044557b4ed9f5ab2c39b84/#fifo8" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo140" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo141" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo144" ; } 
 { { mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo145" ; } 
 { mkfifo "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo146" ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ sort -r -n -m "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo141" "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo145" >"/tmp/pash_1hm9dlJ/275b4d0b7e044557b4ed9f5ab2c39b84/#fifo8" & }
pids_to_kill="${!} ${pids_to_kill}"
{ head -n 1 <"/tmp/pash_1hm9dlJ/275b4d0b7e044557b4ed9f5ab2c39b84/#fifo8" >"/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo146" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo140" -o "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo141" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo144" -o "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo145" -I -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib recv*c25044a9-3ce0-4bb2-a866-572c8b139cc7*1*0*/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo140 recv*f4b19298-04c5-480d-ad97-02a57cefe6f4*1*0*/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo144 & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3.9 aws/s3-put-object.py max-temp/outputs/max-temp-process.sh:temperatures.2015.txt:8:hybrid:max.stdout.txt "/tmp/pash_1hm9dlJ/57426248154648728ee8a0366562834d/#fifo146" $1 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
