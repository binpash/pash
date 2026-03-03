#!/bin/bash
cd $PASH_TOP

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
mkdir -p /tmp/pash_2f9Dvx0/ 
mkdir -p /tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/ 
mkdir -p /tmp/pash_2f9Dvx0/0618050848cd4a33b21a8d7223478225/ 
rm_pash_fifos() {
{ rm -f "/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo15" ; } 
 { { rm -f "/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo16" ; } 
 { { rm -f "/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo21" ; } 
 { { rm -f "/tmp/pash_2f9Dvx0/0618050848cd4a33b21a8d7223478225/#fifo22" ; } 
 { rm -f "/tmp/pash_2f9Dvx0/0618050848cd4a33b21a8d7223478225/#fifo23" ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo15" ; } 
 { { mkfifo "/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo16" ; } 
 { { mkfifo "/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo21" ; } 
 { { mkfifo "/tmp/pash_2f9Dvx0/0618050848cd4a33b21a8d7223478225/#fifo22" ; } 
 { mkfifo "/tmp/pash_2f9Dvx0/0618050848cd4a33b21a8d7223478225/#fifo23" ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_merge "/tmp/pash_2f9Dvx0/0618050848cd4a33b21a8d7223478225/#fifo22" "/tmp/pash_2f9Dvx0/0618050848cd4a33b21a8d7223478225/#fifo23" >"/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo16" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/pashlib-ft recv*346f4c00-ca62-4e46-ba0f-1f398f6eaf43*1*0*/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo15 recv*30d05e46-6454-4f06-be4b-ec5b9da3e97f*1*0*/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo21 & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo15" -o "/tmp/pash_2f9Dvx0/0618050848cd4a33b21a8d7223478225/#fifo22" -I -f -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo21" -o "/tmp/pash_2f9Dvx0/0618050848cd4a33b21a8d7223478225/#fifo23" -I -f -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3 aws/s3-put-object.py ft/nfa-regex-1M.txt "/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo16" $1 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )