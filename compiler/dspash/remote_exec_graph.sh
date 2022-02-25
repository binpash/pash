file=$1

function get_port(){
    comm -23 <(seq 49152 65535 | sort) <(ss -Htan | awk '{print $4}' | cut -d':' -f2 | sort -u) | shuf | head -n 1
}

# pash_redir_output echo "Sending msg to worker manager: $message"
response=($(echo "Exec-Graph: $file" | nc -U "$DSPASH_SOCKET"))
# pash_redir_output echo "Got response from worker manager: $response"   

status=${response[0]} #do something if false
final_output_port=${response[1]}

$PASH_TOP/compiler/dspash/remote_read.sh -l -p $final_output_port