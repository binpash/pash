mkfifo s1 s2

## This way of fixing the problem suffers from some issues.
##
## - First of all, gathering the children after the end of the graph
##   seems to gather more than just the alive nodes. This could lead
##   to killing some random pid in the system. This could potentially
##   be solved by gathering all pids incrementally.
##
## - In addition, this way of getting the last pid does not work if
##   there is more than one output. (This is never the case in our
##   tests, but could be.
##
## - Finally, it is not local, since all of the monitoring happens
##   globally. Ideally, it should be done by a wrapper in each -
##   node. The wrapper should monitor if the node dies, and if so it -
##   should send SIGPIPE to all its producers.

cat ../evaluation/scripts/input/1M.txt > s1 &
echo "Current node: $!"
cat ../evaluation/scripts/input/1M.txt > s2 &
echo "Current node: $!"
cat s1 s2 | head -n 1 &

last=$!

echo "Children pids"
ps --ppid $$ | awk '{print $1}' | grep -E '[0-9]'

echo "Alternative children pids"
jobs -l | awk '{print $1}'

wait $last

echo "Last pid: $last"

ps --ppid $$ | awk '{print $1}' | grep -E '[0-9]' | xargs -n 1 kill -SIGPIPE

rm s1 s2
