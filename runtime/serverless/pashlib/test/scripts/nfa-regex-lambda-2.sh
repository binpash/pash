export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export LOCPATH=/var/task/runtime/locale
export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_2f9Dvx0/ 
mkdir -p /tmp/pash_2f9Dvx0/a65bfe12e9c042c0a483eb706f249675/ 
mkdir -p /tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/ 
rm_pash_fifos() {
{ rm -f "/tmp/pash_2f9Dvx0/a65bfe12e9c042c0a483eb706f249675/#fifo10" ; } 
 { { rm -f "/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo19" ; } 
 { { rm -f "/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo24" ; } 
 { rm -f "/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo25" ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_2f9Dvx0/a65bfe12e9c042c0a483eb706f249675/#fifo10" ; } 
 { { mkfifo "/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo19" ; } 
 { { mkfifo "/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo24" ; } 
 { mkfifo "/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo25" ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_wrap bash -c ' tr A-Z a-z ' <"/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo25" >"/tmp/pash_2f9Dvx0/a65bfe12e9c042c0a483eb706f249675/#fifo10" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_wrap bash -c ' grep "\\(.\\).*\\1\\(.\\).*\\2\\(.\\).*\\3\\(.\\).*\\4" ' <"/tmp/pash_2f9Dvx0/a65bfe12e9c042c0a483eb706f249675/#fifo10" >"/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo19" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo24" -o "/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo25" -I -m 1G -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3 aws/s3-chunk-reader-approx-correction.py "oneliners/inputs/1M.txt" "/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo24" '[{"start": 32768, "end": 65535, "block_id": 1, "shard_id": 0}, {"start": 98304, "end": 131071, "block_id": 3, "shard_id": 1}, {"start": 163840, "end": 196607, "block_id": 5, "shard_id": 2}, {"start": 229376, "end": 262143, "block_id": 7, "shard_id": 3}, {"start": 294912, "end": 327679, "block_id": 9, "shard_id": 4}, {"start": 360448, "end": 393215, "block_id": 11, "shard_id": 5}, {"start": 425984, "end": 458751, "block_id": 13, "shard_id": 6}, {"start": 491520, "end": 524287, "block_id": 15, "shard_id": 7}, {"start": 557056, "end": 589823, "block_id": 17, "shard_id": 8}, {"start": 622592, "end": 655359, "block_id": 19, "shard_id": 9}, {"start": 688128, "end": 720895, "block_id": 21, "shard_id": 10}, {"start": 753664, "end": 786431, "block_id": 23, "shard_id": 11}, {"start": 819200, "end": 851967, "block_id": 25, "shard_id": 12}, {"start": 884736, "end": 917503, "block_id": 27, "shard_id": 13}, {"start": 950272, "end": 983039, "block_id": 29, "shard_id": 14}, {"start": 1015808, "end": 1048576, "block_id": 31, "shard_id": 15}]' shard=1 num_shards=2 job_uid=b4f347ce-735b-4d83-b2f1-2f72e3a32cd7 debug=True window_size=None chunks_per_lambda=16 write_headers=true >"/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo24" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/pashlib-ft send*30d05e46-6454-4f06-be4b-ec5b9da3e97f*0*1*/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo19 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
