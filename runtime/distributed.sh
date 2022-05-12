#!/usr/bin/env bash

# Schematic for distributing computations
# To try server client:

./client.js local 'ls'
./client.js local | jq .
./client.js local 'sleep 5; echo "yo"'
./client.js local | jq .
sleep 5
./client.js local | jq .

# The format is as follows:
# cat FIFOs > $OUT &
#   nc -l -p STAR_RESULT_PORT > FIFO &
#   ./client.js WORKER 'nc -l WORKER_DATA_PORT | PROGRAM | nc -C 158.130.4.212 STAR_RESULT_PORT'
#   cat $IN | nc -N WORKER WORKER_DATA_PORT


cat fifo5555 fifo5556 > "$OUT"
nc -l -p 5555 > fifo5555 &
nc -l -p 5556 > fifo5556 &
# beta runs: `nc -l 5000 | grep -v "onetwo" | tr '[:lower:]' '[:upper:]' | nc -C 158.130.4.212 5555`
./client.js beta 'nc -l 5000 | grep -v "onetwo" | tr "[:lower:]" "[:upper:]" | nc -C 158.130.4.212 5555'
# gamma runs: `nc -l 5000 | grep -v "onetwo" | tr '[:lower:]' '[:upper:]' | nc -C 158.130.4.212 5555`
./client.js gamma 'nc -l 5000 | grep -v "onetwo" | tr "[:lower:]" "[:upper:]" | nc -C 158.130.4.212 5556'
cat "$IN" | nc -N beta.ndr.md 5000
cat "$IN" | nc -N gamma.ndr.md 5000


# Collect results
# Implement: `socat` can listen for multiple connections
nc -l -p 5555 > s1 &
nc -l -p 5556 > s2 &
nc -l -p 5557 > r3 &
nc -l -p 5558 > r4 &

# Things are complicated by the fact that machines default to different
# versions of the BSD netcat (not sure why Debian and Ubuntu default to
# the BSD version of `nc`)
# -N stops after EOF

# Receiver should run
nc -l -p 5000 | tr '[:lower:]' '[:upper:]'  | nc "$DSTAR" 5555

# Distribute load
cat ./a/b | tr 'x' 'x' | nc "$B" 5000
cat ./a/b | tr 'x' 'x' | nc "$C" 5000
cat ./a/b | tr 'x' 'x' | nc "$D" 5000
