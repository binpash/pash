export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export LOCPATH=/var/task/runtime/locale
export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export RUST_BACKTRACE=1
export PASH_ENABLE_RESUMABILITY=1
version=$2
mkdir -p /tmp/pash_mExGzpm/ 
mkdir -p /tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/ 
rm_pash_fifos() {
{ rm -f "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo15" ; } 
 { { rm -f "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo20" ; } 
 { rm -f "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo21" ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo15" ; } 
 { { mkfifo "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo20" ; } 
 { mkfifo "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo21" ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ cat "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo20" | tr A-Z a-z | grep '\(.\).*\1\(.\).*\2' | sort >"/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo15" & }
pids_to_kill="${!} ${pids_to_kill}"
# { runtime/dgsh-tee -i "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo20" -o "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo21" -I -m 1G -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3 aws/s3-chunk-reader-approx-correction.py "oneliners/inputs/100M.txt" "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo20" '[{"start": 3276803, "end": 6553605, "block_id": 1, "shard_id": 0}, {"start": 9830409, "end": 13107211, "block_id": 3, "shard_id": 1}, {"start": 16384015, "end": 19660817, "block_id": 5, "shard_id": 2}, {"start": 22937621, "end": 26214423, "block_id": 7, "shard_id": 3}, {"start": 29491227, "end": 32768029, "block_id": 9, "shard_id": 4}, {"start": 36044833, "end": 39321635, "block_id": 11, "shard_id": 5}, {"start": 42598439, "end": 45875241, "block_id": 13, "shard_id": 6}, {"start": 49152045, "end": 52428847, "block_id": 15, "shard_id": 7}, {"start": 55705651, "end": 58982453, "block_id": 17, "shard_id": 8}, {"start": 62259257, "end": 65536059, "block_id": 19, "shard_id": 9}, {"start": 68812863, "end": 72089665, "block_id": 21, "shard_id": 10}, {"start": 75366469, "end": 78643271, "block_id": 23, "shard_id": 11}, {"start": 81920075, "end": 85196877, "block_id": 25, "shard_id": 12}, {"start": 88473681, "end": 91750483, "block_id": 27, "shard_id": 13}, {"start": 95027287, "end": 98304089, "block_id": 29, "shard_id": 14}, {"start": 101580893, "end": 104857699, "block_id": 31, "shard_id": 15}]' shard=1 num_shards=2 job_uid=74e4fe8f-f4ab-4c17-99a8-6280c4889666 debug=True window_size=None chunks_per_lambda=16 write_headers=false >"/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo20" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/pashlib-ft send*53880cce-bf65-4a27-b383-c4ad58304e7c*0*1*/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo15 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )