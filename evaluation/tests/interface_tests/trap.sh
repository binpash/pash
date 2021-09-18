rm -f trap.txt

trap 'echo fudge > trap.txt' INT

sleep 10
