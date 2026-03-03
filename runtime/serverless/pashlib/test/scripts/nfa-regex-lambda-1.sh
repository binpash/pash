export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export LOCPATH=/var/task/runtime/locale
export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_D06xQF6/ 
mkdir -p /tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/ 
mkdir -p /tmp/pash_D06xQF6/4e97d805a67845bb9d68e34dceeb8543/ 
rm_pash_fifos() {
{ rm -f "/tmp/pash_D06xQF6/4e97d805a67845bb9d68e34dceeb8543/#fifo9" ; } 
 { { rm -f "/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo13" ; } 
 { { rm -f "/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo22" ; } 
 { rm -f "/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo23" ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_D06xQF6/4e97d805a67845bb9d68e34dceeb8543/#fifo9" ; } 
 { { mkfifo "/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo13" ; } 
 { { mkfifo "/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo22" ; } 
 { mkfifo "/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo23" ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_wrap bash -c ' tr A-Z a-z ' <"/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo23" >"/tmp/pash_D06xQF6/4e97d805a67845bb9d68e34dceeb8543/#fifo9" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_wrap bash -c ' grep "\\(.\\).*\\1\\(.\\).*\\2\\(.\\).*\\3\\(.\\).*\\4" ' <"/tmp/pash_D06xQF6/4e97d805a67845bb9d68e34dceeb8543/#fifo9" >"/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo13" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo22" -o "/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo23" -I -m 1G -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3 aws/s3-chunk-reader-approx-correction.py "oneliners/inputs/1M.txt" "/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo22" '[{"start": 0, "end": 524287, "block_id": 0, "shard_id": 0}]' shard=0 num_shards=2 job_uid=60076533-7e02-4ff9-aca9-bc9a06a4a7b0 debug=True window_size=None chunks_per_lambda=1 write_headers=true >"/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo22" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/pashlib-ft send*0043f988-7b29-4455-abee-30c80ae67b98*0*1*/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo13 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
