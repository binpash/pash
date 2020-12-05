#!/bin/bash

## Kill the original script if the trap is caught
cleanup()
{
    echo seq_pid: $seq_pid ;
    kill -9 $seq_pid ;
    exit 0
}

## Create a trap that handles SIGUSR1 and kills the original script execution
echo "Trap set in $BASHPID"
trap cleanup SIGUSR1

## Execute the original script
echo $BASHPID
./exp2.sh &
seq_pid=$!
wait

echo "popo"