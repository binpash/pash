#!/bin/bash

## Two possible implementations:
## 1. Either create a trap and when compilation ends signal the original shell to kill the 
##    sequential and execute the parallel
## 2. Either wait for both to end and if the sequential ends first then just kill the compilation, 
##    else just kill the sequential and run the parallel.

## Kill the original script if the trap is caught
cleanup()
{
    echo seq_pid: $seq_pid ;
    kill -9 $seq_pid ;
    exit 0
}

still_alive()
{
    jobs -p | tr '\n' ' '
}

## Create a trap that handles SIGUSR1 and kills the original script execution
echo "Trap set in $BASHPID"
trap cleanup SIGUSR1

## Execute the original script
echo $BASHPID
./sleep.sh 10 &
seq_pid=$!

./sleep.sh 1 &
comp_pid=$!

## What if they both terminate before wait
echo "Still alive: $(still_alive)"
wait -n $seq_pid $comp_pid
alive_pid=$(still_alive)
echo "Still alive: $(still_alive)"
## If the sequential is still alive we want to kill it
if [ "$seq_pid" -eq "$alive_pid" ]; then
    echo "Compilation was done first!"
    kill -n 9 "$seq_pid" 2> /dev/null
    kill_status=$?
    ## If kill failed it means it was already dead
    if [ "$kill_status" -eq 0 ]; then
        ## TODO: Run the parallel script
        echo "Run parallel"
    fi
else
    echo "Sequential was done first!"
    ## If this fails (meaning that compilation is done) we do not care
    kill -n 9 "$comp_pid" 2> /dev/null
fi

echo "Done"