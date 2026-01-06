export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_vtZsjBm/ 
mkdir -p /tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/ 
mkdir -p /tmp/pash_vtZsjBm/68a7f4eeec164679bf1bc1313be933aa/ 
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
{ rm -f "/tmp/pash_vtZsjBm/68a7f4eeec164679bf1bc1313be933aa/#fifo2" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo129" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo133" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo137" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo141" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo145" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo149" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo153" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo157" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo161" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo165" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo169" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo173" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo177" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo181" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo185" ; } 
 { { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo189" ; } 
 { rm -f "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo388" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_vtZsjBm/68a7f4eeec164679bf1bc1313be933aa/#fifo2" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo129" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo133" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo137" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo141" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo145" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo149" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo153" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo157" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo161" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo165" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo169" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo173" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo177" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo181" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo185" ; } 
 { { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo189" ; } 
 { mkfifo "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo388" ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ cat "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo388" >"/tmp/pash_vtZsjBm/68a7f4eeec164679bf1bc1313be933aa/#fifo2" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_split "/tmp/pash_vtZsjBm/68a7f4eeec164679bf1bc1313be933aa/#fifo2" 1000000 "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo129" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo133" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo137" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo141" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo145" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo149" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo153" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo157" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo161" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo165" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo169" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo173" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo177" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo181" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo185" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo189" & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3.9 aws/s3-get-object.py "covid-mts/inputs/in.csv" "/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo388" & }
pids_to_kill="${!} ${pids_to_kill}"
{ /opt/pashlib send*1a6780dc-a4ec-4674-ada3-bb15a795f216*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo129 send*20d35c51-b875-4f3d-b557-095ac6e111c0*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo133 send*20ab305d-7440-49e9-8fdb-8c81a9c65180*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo137 send*f81db0cf-92c5-4293-b4a0-17efc6263e3b*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo141 send*087f25fc-0dc3-4ea3-90fa-80629754a15c*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo145 send*65552f6e-bb3d-4b72-893f-b301e1d957c3*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo149 send*79a63115-19a8-4464-9242-3d3dd452a5e3*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo153 send*e4d37939-2bab-4e28-8e9b-7affb0c7bce4*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo157 send*68d89b08-deb3-440e-b83c-748ec35d0db1*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo161 send*09355f6f-191f-43df-9c9a-235e7ce3df05*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo165 send*2ddcda7d-3d29-43e1-bd87-c391bb157dfb*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo169 send*f5996726-0dbf-4cd1-b1f6-8c58eb374902*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo173 send*6e27162f-a952-41aa-9fc1-0b2392979fdc*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo177 send*c0ea277c-95ef-499a-a1a8-af0dc46b11b0*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo181 send*7b108f7b-a6e5-4504-8a27-56d21dbbc868*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo185 send*b6e5a89e-789e-4f99-8b96-2fe338ed6859*0*1*/tmp/pash_vtZsjBm/5fe2a53b92b64096a641533fb280941f/#fifo189 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
