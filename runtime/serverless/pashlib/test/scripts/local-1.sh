#!/bin/bash
cd $PASH_TOP

# Kill entire process group on exit
cleanup() {
    trap - EXIT INT TERM
    kill -- -$$ 2>/dev/null || true
}
trap cleanup EXIT INT TERM

export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export LOCPATH=/var/task/runtime/locale
export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_mExGzpm/ 
mkdir -p /tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/ 
rm_pash_fifos() {
{ rm -f "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo9" ; } 
 { { rm -f "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo18" ; } 
 { rm -f "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo19" ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo9" ; } 
 { { mkfifo "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo18" ; } 
 { mkfifo "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo19" ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ sort <"/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo19" >"/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo9" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo18" -o "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo19" -I -m 1G -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3 aws/s3-get-object.py "oneliners/inputs/1M.txt" "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo18" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/pashlib-ft send*12329258-0e5c-4b3f-a37e-613d187b9961*0*1*/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo9 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )