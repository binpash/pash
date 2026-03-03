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
{ rm -f "/tmp/pash_D06xQF6/4e97d805a67845bb9d68e34dceeb8543/#fifo10" ; } 
 { { rm -f "/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo19" ; } 
 { { rm -f "/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo24" ; } 
 { rm -f "/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo25" ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_D06xQF6/4e97d805a67845bb9d68e34dceeb8543/#fifo10" ; } 
 { { mkfifo "/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo19" ; } 
 { { mkfifo "/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo24" ; } 
 { mkfifo "/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo25" ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_wrap bash -c ' tr A-Z a-z ' <"/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo25" >"/tmp/pash_D06xQF6/4e97d805a67845bb9d68e34dceeb8543/#fifo10" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_wrap bash -c ' grep "\\(.\\).*\\1\\(.\\).*\\2\\(.\\).*\\3\\(.\\).*\\4" ' <"/tmp/pash_D06xQF6/4e97d805a67845bb9d68e34dceeb8543/#fifo10" >"/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo19" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo24" -o "/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo25" -I -m 1G -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3 aws/s3-chunk-reader-approx-correction.py "oneliners/inputs/1M.txt" "/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo24" '[{"start": 524288, "end": 1048576, "block_id": 1, "shard_id": 0}]' shard=1 num_shards=2 job_uid=60076533-7e02-4ff9-aca9-bc9a06a4a7b0 debug=True window_size=None chunks_per_lambda=1 write_headers=true >"/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo24" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/pashlib-ft send*bf85c55a-06f3-4b14-a90b-79fe17d56830*0*1*/tmp/pash_D06xQF6/015e127d417245998600556616ca0b8f/#fifo19 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
