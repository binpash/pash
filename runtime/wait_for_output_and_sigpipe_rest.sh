#!/usr/bin/env bash

## TODO: Give it the output pid as an argument
wait $@

## It is assumed that $distro is set when this is called.

## Note: We need the || true after the grep so that it doesn't exit with error if it finds nothing.

# now do different things depending on distro
case "$distro" in
    freebsd*)  
        # not sure at all about this one
        pids_to_kill="$(ps -efl $BASHPID |awk '{print $1}' | { grep -E '[0-9]' || true; } )"
        ;;
    *)
        pids_to_kill="$(ps --ppid $BASHPID |awk '{print $1}' | { grep -E '[0-9]' || true; } )"
        ;;
esac

for pid in $pids_to_kill
do
    (> /dev/null 2>&1 kill -SIGPIPE $pid || true)
done

## Only wait if that is not empty
if [ ! -z "$var" ];
then
    ## Need to wait to not leave any zombies
    wait
fi
