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
{ rm -f "/tmp/pash_f7NsyHc/9805ed4004cf44a88d5e3f51f2da8f23/#fifo11" ; } 
 { { rm -f "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo19" ; } 
 { { rm -f "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo34" ; } 
 { rm -f "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo35" ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_f7NsyHc/9805ed4004cf44a88d5e3f51f2da8f23/#fifo11" ; } 
 { { mkfifo "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo19" ; } 
 { { mkfifo "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo34" ; } 
 { mkfifo "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo35" ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_wrap bash -c ' tr A-Z a-z ' <"/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo35" >"/tmp/pash_f7NsyHc/9805ed4004cf44a88d5e3f51f2da8f23/#fifo11" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_wrap bash -c ' grep "\\(.\\).*\\1\\(.\\).*\\2" ' <"/tmp/pash_f7NsyHc/9805ed4004cf44a88d5e3f51f2da8f23/#fifo11" >"/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo19" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo34" -o "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo35" -I -m 1G -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3 aws/s3-chunk-reader-approx-correction.py "oneliners/inputs/100M.txt" "/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo34" '[{"start": 0, "end": 1638400, "block_id": 0, "shard_id": 0}, {"start": 6553604, "end": 8192004, "block_id": 4, "shard_id": 1}, {"start": 13107208, "end": 14745608, "block_id": 8, "shard_id": 2}, {"start": 19660812, "end": 21299212, "block_id": 12, "shard_id": 3}, {"start": 26214416, "end": 27852816, "block_id": 16, "shard_id": 4}, {"start": 32768020, "end": 34406420, "block_id": 20, "shard_id": 5}, {"start": 39321624, "end": 40960024, "block_id": 24, "shard_id": 6}, {"start": 45875228, "end": 47513628, "block_id": 28, "shard_id": 7}, {"start": 52428832, "end": 54067232, "block_id": 32, "shard_id": 8}, {"start": 58982436, "end": 60620836, "block_id": 36, "shard_id": 9}, {"start": 65536040, "end": 67174440, "block_id": 40, "shard_id": 10}, {"start": 72089644, "end": 73728044, "block_id": 44, "shard_id": 11}, {"start": 78643248, "end": 80281648, "block_id": 48, "shard_id": 12}, {"start": 85196852, "end": 86835252, "block_id": 52, "shard_id": 13}, {"start": 91750456, "end": 93388856, "block_id": 56, "shard_id": 14}, {"start": 98304060, "end": 99942460, "block_id": 60, "shard_id": 15}]' shard=0 num_shards=4 job_uid=2a008354-1b6a-438a-908a-7abf1f061794 debug=True window_size=None chunks_per_lambda=16 write_headers=true >"/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo34" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/pashlib-ft send*7746ca28-2248-41b6-8736-cb943021e565*0*1*/tmp/pash_f7NsyHc/4a0e9f979abb43acbcdb621b0b5b67cd/#fifo19 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )