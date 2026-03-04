#!/bin/bash
cd $PASH_TOP
export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export LOCPATH=/var/task/runtime/locale
export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_f7NsyHc/ 
mkdir -p /tmp/pash_f7NsyHc/12ef22836c754348942f07eb684e81cb/ 
mkdir -p /tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/ 

rm_pash_fifos() {
{ rm -f "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo21" ; } 
 { { rm -f "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo22" ; } 
 { { rm -f "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo27" ; } 
 { { rm -f "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo30" ; } 
 { { rm -f "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo33" ; } 
 { { rm -f "/tmp/pash_f7NsyHc/12ef22836c754348942f07eb684e81cb/#fifo34" ; } 
 { { rm -f "/tmp/pash_f7NsyHc/12ef22836c754348942f07eb684e81cb/#fifo35" ; } 
 { { rm -f "/tmp/pash_f7NsyHc/12ef22836c754348942f07eb684e81cb/#fifo36" ; } 
 { rm -f "/tmp/pash_f7NsyHc/12ef22836c754348942f07eb684e81cb/#fifo37" ; } ; } ; } ; } ; } ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo21" ; } 
 { { mkfifo "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo22" ; } 
 { { mkfifo "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo27" ; } 
 { { mkfifo "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo30" ; } 
 { { mkfifo "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo33" ; } 
 { { mkfifo "/tmp/pash_f7NsyHc/12ef22836c754348942f07eb684e81cb/#fifo34" ; } 
 { { mkfifo "/tmp/pash_f7NsyHc/12ef22836c754348942f07eb684e81cb/#fifo35" ; } 
 { { mkfifo "/tmp/pash_f7NsyHc/12ef22836c754348942f07eb684e81cb/#fifo36" ; } 
 { mkfifo "/tmp/pash_f7NsyHc/12ef22836c754348942f07eb684e81cb/#fifo37" ; } ; } ; } ; } ; } ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_merge "/tmp/pash_f7NsyHc/12ef22836c754348942f07eb684e81cb/#fifo34" "/tmp/pash_f7NsyHc/12ef22836c754348942f07eb684e81cb/#fifo35" "/tmp/pash_f7NsyHc/12ef22836c754348942f07eb684e81cb/#fifo36" "/tmp/pash_f7NsyHc/12ef22836c754348942f07eb684e81cb/#fifo37" >"/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo22" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/pashlib-ft recv*7746ca28-2248-41b6-8736-cb943021e565*1*0*/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo21 recv*cdb65968-53c7-4582-8531-e4c83d3002b8*1*0*/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo27 recv*c8a71b12-c941-41c6-8a57-c51c573eee02*1*0*/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo30 recv*92791fda-fba6-4138-b184-e171ba2425f1*1*0*/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo33 & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo21" -o "/tmp/pash_f7NsyHc/12ef22836c754348942f07eb684e81cb/#fifo34" -I -f -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo27" -o "/tmp/pash_f7NsyHc/12ef22836c754348942f07eb684e81cb/#fifo35" -I -f -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo30" -o "/tmp/pash_f7NsyHc/12ef22836c754348942f07eb684e81cb/#fifo36" -I -f -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo33" -o "/tmp/pash_f7NsyHc/12ef22836c754348942f07eb684e81cb/#fifo37" -I -f -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3 aws/s3-put-object.py ft/nfa-regex-100M.txt "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo22" $1 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )