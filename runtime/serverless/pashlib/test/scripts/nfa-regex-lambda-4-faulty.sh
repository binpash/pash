export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export LOCPATH=/var/task/runtime/locale
export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_f7NsyHc/ 
mkdir -p /tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/ 
mkdir -p /tmp/pash_f7NsyHc/9805ed4004cf44a88d5e3f51f2da8f23/ 

rm_pash_fifos() {
{ rm -f "/tmp/pash_f7NsyHc/9805ed4004cf44a88d5e3f51f2da8f23/#fifo14" ; } 
 { { rm -f "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo31" ; } 
 { { rm -f "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo40" ; } 
 { rm -f "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo41" ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_f7NsyHc/9805ed4004cf44a88d5e3f51f2da8f23/#fifo14" ; } 
 { { mkfifo "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo31" ; } 
 { { mkfifo "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo40" ; } 
 { mkfifo "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo41" ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_wrap bash -c ' tr A-Z a-z ' <"/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo41" >"/tmp/pash_f7NsyHc/9805ed4004cf44a88d5e3f51f2da8f23/#fifo14" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_wrap bash -c ' grep "\\(.\\).*\\1\\(.\\).*\\2" ' <"/tmp/pash_f7NsyHc/9805ed4004cf44a88d5e3f51f2da8f23/#fifo14" >"/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo31" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo40" -o "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo41" -I -m 1G -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3 aws/s3-chunk-reader-approx-correction.py "oneliners/inputs/100M.txt" "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo40" '[{"start": 4915203, "end": 6553603, "block_id": 3, "shard_id": 0}, {"start": 11468807, "end": 13107207, "block_id": 7, "shard_id": 1}, {"start": 18022411, "end": 19660811, "block_id": 11, "shard_id": 2}, {"start": 24576015, "end": 26214415, "block_id": 15, "shard_id": 3}, {"start": 31129619, "end": 32768019, "block_id": 19, "shard_id": 4}, {"start": 37683223, "end": 39321623, "block_id": 23, "shard_id": 5}, {"start": 44236827, "end": 45875227, "block_id": 27, "shard_id": 6}, {"start": 50790431, "end": 52428831, "block_id": 31, "shard_id": 7}, {"start": 57344035, "end": 58982435, "block_id": 35, "shard_id": 8}, {"start": 63897639, "end": 65536039, "block_id": 39, "shard_id": 9}, {"start": 70451243, "end": 72089643, "block_id": 43, "shard_id": 10}, {"start": 77004847, "end": 78643247, "block_id": 47, "shard_id": 11}, {"start": 83558451, "end": 85196851, "block_id": 51, "shard_id": 12}, {"start": 90112055, "end": 91750455, "block_id": 55, "shard_id": 13}, {"start": 96665659, "end": 98304059, "block_id": 59, "shard_id": 14}, {"start": 103219263, "end": 104857699, "block_id": 63, "shard_id": 15}]' shard=3 num_shards=4 job_uid=2a008354-1b6a-438a-908a-7abf1f061794 debug=True window_size=None chunks_per_lambda=16 write_headers=true >"/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo40" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/pashlib-ft send*92791fda-fba6-4138-b184-e171ba2425f1*0*1*/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo31 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )