cmd=$1

# function get_ports(){
#     comm -23 <(seq 49152 65535 | sort) <(ss -Htan | awk '{print $4}' | cut -d':' -f2 | sort -u) | shuf | head -n 2
# }

# ports=($(get_ports))
# from_port=${ports[0]}
# to_port=${ports[1]}

# pash_redir_output echo "Sending msg to worker manager: $message"
response=($(echo "Get-Worker: $cmd" | nc -U "$DSPASH_SOCKET"))
# pash_redir_output echo "Got response from worker manager: $response"   

status=${response[0]} #do something if false
ip=${response[1]}
from_port=${response[2]}
to_port=${response[3]}

nc -d $ip $to_port &
nc $ip $from_port