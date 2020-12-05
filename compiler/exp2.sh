# # # !/bin/bash

# trap "echo hi > temp ; exit 0" SIGUSR1


echo "Hi! I am the sleeper: $BASHPID"
sleep 10000
# wait