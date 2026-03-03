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
{ rm -f "/tmp/pash_2f9Dvx0/a65bfe12e9c042c0a483eb706f249675/#fifo9" ; } 
 { { rm -f "/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo13" ; } 
 { { rm -f "/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo22" ; } 
 { rm -f "/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo23" ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_2f9Dvx0/a65bfe12e9c042c0a483eb706f249675/#fifo9" ; } 
 { { mkfifo "/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo13" ; } 
 { { mkfifo "/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo22" ; } 
 { mkfifo "/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo23" ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_wrap bash -c ' tr A-Z a-z ' <"/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo23" >"/tmp/pash_2f9Dvx0/a65bfe12e9c042c0a483eb706f249675/#fifo9" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_wrap bash -c ' grep "\\(.\\).*\\1\\(.\\).*\\2\\(.\\).*\\3\\(.\\).*\\4" ' <"/tmp/pash_2f9Dvx0/a65bfe12e9c042c0a483eb706f249675/#fifo9" >"/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo13" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo22" -o "/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo23" -I -m 1G -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3 aws/s3-chunk-reader-approx-correction.py "oneliners/inputs/1M.txt" "/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo22" '[{"start": 0, "end": 32767, "block_id": 0, "shard_id": 0}, {"start": 65536, "end": 98303, "block_id": 2, "shard_id": 1}, {"start": 131072, "end": 163839, "block_id": 4, "shard_id": 2}, {"start": 196608, "end": 229375, "block_id": 6, "shard_id": 3}, {"start": 262144, "end": 294911, "block_id": 8, "shard_id": 4}, {"start": 327680, "end": 360447, "block_id": 10, "shard_id": 5}, {"start": 393216, "end": 425983, "block_id": 12, "shard_id": 6}, {"start": 458752, "end": 491519, "block_id": 14, "shard_id": 7}, {"start": 524288, "end": 557055, "block_id": 16, "shard_id": 8}, {"start": 589824, "end": 622591, "block_id": 18, "shard_id": 9}, {"start": 655360, "end": 688127, "block_id": 20, "shard_id": 10}, {"start": 720896, "end": 753663, "block_id": 22, "shard_id": 11}, {"start": 786432, "end": 819199, "block_id": 24, "shard_id": 12}, {"start": 851968, "end": 884735, "block_id": 26, "shard_id": 13}, {"start": 917504, "end": 950271, "block_id": 28, "shard_id": 14}, {"start": 983040, "end": 1015807, "block_id": 30, "shard_id": 15}]' shard=0 num_shards=2 job_uid=b4f347ce-735b-4d83-b2f1-2f72e3a32cd7 debug=True window_size=None chunks_per_lambda=16 write_headers=true >"/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo22" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/pashlib-ft send*346f4c00-ca62-4e46-ba0f-1f398f6eaf43*0*1*/tmp/pash_2f9Dvx0/0e8e168fe76442fa9a8938c7fa76e431/#fifo13 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )