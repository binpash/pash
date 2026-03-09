
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
{ rm -f "/tmp/pash_f7NsyHc/9805ed4004cf44a88d5e3f51f2da8f23/#fifo12" ; } 
 { { rm -f "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo25" ; } 
 { { rm -f "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo36" ; } 
 { rm -f "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo37" ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_f7NsyHc/9805ed4004cf44a88d5e3f51f2da8f23/#fifo12" ; } 
 { { mkfifo "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo25" ; } 
 { { mkfifo "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo36" ; } 
 { mkfifo "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo37" ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_wrap bash -c ' tr A-Z a-z ' <"/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo37" >"/tmp/pash_f7NsyHc/9805ed4004cf44a88d5e3f51f2da8f23/#fifo12" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_wrap bash -c ' grep "\\(.\\).*\\1\\(.\\).*\\2" ' <"/tmp/pash_f7NsyHc/9805ed4004cf44a88d5e3f51f2da8f23/#fifo12" >"/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo25" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo36" -o "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo37" -I -m 1G -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3 aws/s3-chunk-reader-approx-correction.py "oneliners/inputs/100M.txt" "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo36" '[{"start": 1638401, "end": 3276801, "block_id": 1, "shard_id": 0}, {"start": 8192005, "end": 9830405, "block_id": 5, "shard_id": 1}, {"start": 14745609, "end": 16384009, "block_id": 9, "shard_id": 2}, {"start": 21299213, "end": 22937613, "block_id": 13, "shard_id": 3}, {"start": 27852817, "end": 29491217, "block_id": 17, "shard_id": 4}, {"start": 34406421, "end": 36044821, "block_id": 21, "shard_id": 5}, {"start": 40960025, "end": 42598425, "block_id": 25, "shard_id": 6}, {"start": 47513629, "end": 49152029, "block_id": 29, "shard_id": 7}, {"start": 54067233, "end": 55705633, "block_id": 33, "shard_id": 8}, {"start": 60620837, "end": 62259237, "block_id": 37, "shard_id": 9}, {"start": 67174441, "end": 68812841, "block_id": 41, "shard_id": 10}, {"start": 73728045, "end": 75366445, "block_id": 45, "shard_id": 11}, {"start": 80281649, "end": 81920049, "block_id": 49, "shard_id": 12}, {"start": 86835253, "end": 88473653, "block_id": 53, "shard_id": 13}, {"start": 93388857, "end": 95027257, "block_id": 57, "shard_id": 14}, {"start": 99942461, "end": 101580861, "block_id": 61, "shard_id": 15}]' shard=1 num_shards=4 job_uid=2a008354-1b6a-438a-908a-7abf1f061794 debug=True window_size=None chunks_per_lambda=16 write_headers=true >"/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo36" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/pashlib-ft send*cdb65968-53c7-4582-8531-e4c83d3002b8*0*1*/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo25 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )