Host=$1
fromPort=$2
toPort=$3
cmd=$4

### failed to use one port for both reading and writing.
# The stream doesn't seem to end even with -d
# mkfifo $readfifo $writefifo
# cat $writefifo | nc $Host 5000 > $readfifo &
# cat $readfifo | $cmd > $writefifo
# rm $readfifo $writefifo
nc -N -l -d $Host $fromPort | bash -c "$cmd" | nc -l -N $Host $toPort
# nc -d $Host $fromPort | bash -c "$cmd" | nc $Host $toPort