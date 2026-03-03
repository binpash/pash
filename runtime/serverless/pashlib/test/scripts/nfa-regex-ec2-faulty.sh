#!/bin/bash
cd $PASH_TOP
export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export LOCPATH=/var/task/runtime/locale
export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_Rx37Zq5/ 
mkdir -p /tmp/pash_Rx37Zq5/2b1ce94b3fd34c70ba7a89753dd3667d/ 
mkdir -p /tmp/pash_Rx37Zq5/417af50694ad42819b8311edb72c13c3/ 
rm_pash_fifos() {
{ rm -f "/tmp/pash_Rx37Zq5/417af50694ad42819b8311edb72c13c3/#fifo15" ; } 
 { { rm -f "/tmp/pash_Rx37Zq5/417af50694ad42819b8311edb72c13c3/#fifo16" ; } 
 { { rm -f "/tmp/pash_Rx37Zq5/417af50694ad42819b8311edb72c13c3/#fifo21" ; } 
 { { rm -f "/tmp/pash_Rx37Zq5/2b1ce94b3fd34c70ba7a89753dd3667d/#fifo22" ; } 
 { rm -f "/tmp/pash_Rx37Zq5/2b1ce94b3fd34c70ba7a89753dd3667d/#fifo23" ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_Rx37Zq5/417af50694ad42819b8311edb72c13c3/#fifo15" ; } 
 { { mkfifo "/tmp/pash_Rx37Zq5/417af50694ad42819b8311edb72c13c3/#fifo16" ; } 
 { { mkfifo "/tmp/pash_Rx37Zq5/417af50694ad42819b8311edb72c13c3/#fifo21" ; } 
 { { mkfifo "/tmp/pash_Rx37Zq5/2b1ce94b3fd34c70ba7a89753dd3667d/#fifo22" ; } 
 { mkfifo "/tmp/pash_Rx37Zq5/2b1ce94b3fd34c70ba7a89753dd3667d/#fifo23" ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_merge "/tmp/pash_Rx37Zq5/2b1ce94b3fd34c70ba7a89753dd3667d/#fifo22" "/tmp/pash_Rx37Zq5/2b1ce94b3fd34c70ba7a89753dd3667d/#fifo23" >"/tmp/pash_Rx37Zq5/417af50694ad42819b8311edb72c13c3/#fifo16" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/pashlib-ft recv*4fa45390-6e5b-4123-963e-786edd2f8586*1*0*/tmp/pash_Rx37Zq5/417af50694ad42819b8311edb72c13c3/#fifo15 recv*a69a059b-dc9d-434b-bc2a-192eb626a379*1*0*/tmp/pash_Rx37Zq5/417af50694ad42819b8311edb72c13c3/#fifo21 & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_Rx37Zq5/417af50694ad42819b8311edb72c13c3/#fifo15" -o "/tmp/pash_Rx37Zq5/2b1ce94b3fd34c70ba7a89753dd3667d/#fifo22" -I -f -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_Rx37Zq5/417af50694ad42819b8311edb72c13c3/#fifo21" -o "/tmp/pash_Rx37Zq5/2b1ce94b3fd34c70ba7a89753dd3667d/#fifo23" -I -f -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3 aws/s3-put-object.py ft/nfa-regex-100M.txt "/tmp/pash_Rx37Zq5/417af50694ad42819b8311edb72c13c3/#fifo16" $1 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )