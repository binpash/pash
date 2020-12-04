#!/bin/bash

echo "Trap set in $BASHPID"
trap "echo hi > temp ; exit 0" SIGUSR1

echo $BASHPID
sleep 10000 &
wait
