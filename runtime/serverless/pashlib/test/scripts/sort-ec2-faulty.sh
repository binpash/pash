#!/bin/bash
cd $PASH_TOP

export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export LOCPATH=/var/task/runtime/locale
export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_NSonVV8/ 
mkdir -p /tmp/pash_NSonVV8/706abff955484bd5a54470191dee43d0/ 
mkdir -p /tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/ 
rm_pash_fifos() {
{ rm -f "/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo11" ; } 
 { { rm -f "/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo12" ; } 
 { { rm -f "/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo17" ; } 
 { { rm -f "/tmp/pash_NSonVV8/706abff955484bd5a54470191dee43d0/#fifo18" ; } 
 { rm -f "/tmp/pash_NSonVV8/706abff955484bd5a54470191dee43d0/#fifo19" ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo11" ; } 
 { { mkfifo "/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo12" ; } 
 { { mkfifo "/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo17" ; } 
 { { mkfifo "/tmp/pash_NSonVV8/706abff955484bd5a54470191dee43d0/#fifo18" ; } 
 { mkfifo "/tmp/pash_NSonVV8/706abff955484bd5a54470191dee43d0/#fifo19" ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ cat "/tmp/pash_NSonVV8/706abff955484bd5a54470191dee43d0/#fifo18" "/tmp/pash_NSonVV8/706abff955484bd5a54470191dee43d0/#fifo19" >"/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo12" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/pashlib-ft recv*4c1a5c79-5201-459a-83ae-1a8225370842*1*0*/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo11 recv*b04e1ae6-5eb1-43b1-95a0-7a8334174ba1*1*0*/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo17 & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo11" -o "/tmp/pash_NSonVV8/706abff955484bd5a54470191dee43d0/#fifo18" -I -f -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo17" -o "/tmp/pash_NSonVV8/706abff955484bd5a54470191dee43d0/#fifo19" -I -f -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3 aws/s3-put-object.py ft/stateful-faulty.txt "/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo12" $1 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
