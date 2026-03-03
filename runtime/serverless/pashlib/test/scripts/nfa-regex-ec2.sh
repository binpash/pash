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
mkdir -p /tmp/pash_D06xQF6/ 
mkdir -p /tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/ 
mkdir -p /tmp/pash_D06xQF6/79573b5630e245b38569737e0c39c50b/ 
rm_pash_fifos() {
{ rm -f "/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo15" ; } 
 { { rm -f "/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo16" ; } 
 { { rm -f "/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo21" ; } 
 { { rm -f "/tmp/pash_D06xQF6/79573b5630e245b38569737e0c39c50b/#fifo22" ; } 
 { rm -f "/tmp/pash_D06xQF6/79573b5630e245b38569737e0c39c50b/#fifo23" ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo15" ; } 
 { { mkfifo "/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo16" ; } 
 { { mkfifo "/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo21" ; } 
 { { mkfifo "/tmp/pash_D06xQF6/79573b5630e245b38569737e0c39c50b/#fifo22" ; } 
 { mkfifo "/tmp/pash_D06xQF6/79573b5630e245b38569737e0c39c50b/#fifo23" ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_merge "/tmp/pash_D06xQF6/79573b5630e245b38569737e0c39c50b/#fifo22" "/tmp/pash_D06xQF6/79573b5630e245b38569737e0c39c50b/#fifo23" >"/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo16" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/pashlib-ft recv*0043f988-7b29-4455-abee-30c80ae67b98*1*0*/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo15 recv*bf85c55a-06f3-4b14-a90b-79fe17d56830*1*0*/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo21 & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo15" -o "/tmp/pash_D06xQF6/79573b5630e245b38569737e0c39c50b/#fifo22" -I -f -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo21" -o "/tmp/pash_D06xQF6/79573b5630e245b38569737e0c39c50b/#fifo23" -I -f -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3 aws/s3-put-object.py ft/nfa-regex-1M.txt "/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo16" $1 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )