ir_file=$1

# pash_redir_output echo "Sending msg to worker manager: $message"
response=($(echo "Exec-Graph: $ir_file $declared_functions" | nc -U "$DSPASH_SOCKET"))
# pash_redir_output echo "Got response from worker manager: $response"   

status=${response[0]} #do something if false
script_to_execute=${response[1]}

source "$script_to_execute"
