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
{ rm -f "/tmp/pash_f7NsyHc/9805ed4004cf44a88d5e3f51f2da8f23/#fifo13" ; } 
 { { rm -f "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo28" ; } 
 { { rm -f "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo38" ; } 
 { rm -f "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo39" ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_f7NsyHc/9805ed4004cf44a88d5e3f51f2da8f23/#fifo13" ; } 
 { { mkfifo "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo28" ; } 
 { { mkfifo "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo38" ; } 
 { mkfifo "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo39" ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_wrap bash -c ' tr A-Z a-z ' <"/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo39" >"/tmp/pash_f7NsyHc/9805ed4004cf44a88d5e3f51f2da8f23/#fifo13" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_wrap bash -c ' grep "\\(.\\).*\\1\\(.\\).*\\2" ' <"/tmp/pash_f7NsyHc/9805ed4004cf44a88d5e3f51f2da8f23/#fifo13" >"/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo28" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo38" -o "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo39" -I -m 1G -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3 aws/s3-chunk-reader-approx-correction.py "oneliners/inputs/100M.txt" "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo38" '[{"start": 3276802, "end": 4915202, "block_id": 2, "shard_id": 0}, {"start": 9830406, "end": 11468806, "block_id": 6, "shard_id": 1}, {"start": 16384010, "end": 18022410, "block_id": 10, "shard_id": 2}, {"start": 22937614, "end": 24576014, "block_id": 14, "shard_id": 3}, {"start": 29491218, "end": 31129618, "block_id": 18, "shard_id": 4}, {"start": 36044822, "end": 37683222, "block_id": 22, "shard_id": 5}, {"start": 42598426, "end": 44236826, "block_id": 26, "shard_id": 6}, {"start": 49152030, "end": 50790430, "block_id": 30, "shard_id": 7}, {"start": 55705634, "end": 57344034, "block_id": 34, "shard_id": 8}, {"start": 62259238, "end": 63897638, "block_id": 38, "shard_id": 9}, {"start": 68812842, "end": 70451242, "block_id": 42, "shard_id": 10}, {"start": 75366446, "end": 77004846, "block_id": 46, "shard_id": 11}, {"start": 81920050, "end": 83558450, "block_id": 50, "shard_id": 12}, {"start": 88473654, "end": 90112054, "block_id": 54, "shard_id": 13}, {"start": 95027258, "end": 96665658, "block_id": 58, "shard_id": 14}, {"start": 101580862, "end": 103219262, "block_id": 62, "shard_id": 15}]' shard=2 num_shards=4 job_uid=2a008354-1b6a-438a-908a-7abf1f061794 debug=True window_size=None chunks_per_lambda=16 write_headers=true >"/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo38" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/pashlib-ft send*c8a71b12-c941-41c6-8a57-c51c573eee02*0*1*/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo28 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )