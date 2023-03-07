#!/bin/bash

## Contains utility function that are used by multiple pash frontends

pash_wait_until_unix_socket_listening()
{
    local server_name=$1
    local socket=$2
    ## Only wait for a limited amount of time.
    ## If the daemon cannot start listening in ~ 1 second,
    ##   then it must have crashed or so.
    i=0
    ## This is a magic number to make sure that we wait enough
    maximum_retries=1000
    ## For some reason, `nc -z` doesn't work on livestar (it always returns error)
    ## and therefore we need to send something. 
    until  echo "Daemon Start" 2> /dev/null | nc -U "$socket" >/dev/null 2>&1 ; 
    do 
        ## TODO: Can we wait for the daemon in a better way?
        sleep 0.01
        i=$((i+1))
        if [ $i -eq $maximum_retries ]; then
            echo "Error: Maximum retries: $maximum_retries exceeded when waiting for server: ${server_name} to bind to socket: ${socket}!" 1>&2
            echo "Exiting..." 1>&2
            exit 1
        fi
    done
}

pash_communicate_unix_socket()
{
    local server_name=$1
    local socket=$2
    local message=$3
    pash_redir_output echo "Sending msg to ${server_name}: $message"
    daemon_response=$(echo "$message" | nc -U "${socket}")
    pash_redir_output echo "Got response from ${server_name}: $daemon_response"
    echo "$daemon_response"
}

export -f pash_wait_until_unix_socket_listening
export -f pash_communicate_unix_socket
