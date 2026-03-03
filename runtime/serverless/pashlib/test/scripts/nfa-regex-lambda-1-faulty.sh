export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export LOCPATH=/var/task/runtime/locale
export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_Rx37Zq5/ 
mkdir -p /tmp/pash_Rx37Zq5/9a07e1830b874c71989ee6e88728c8ff/ 
mkdir -p /tmp/pash_Rx37Zq5/417af50694ad42819b8311edb72c13c3/ 
rm_pash_fifos() {
{ rm -f "/tmp/pash_Rx37Zq5/9a07e1830b874c71989ee6e88728c8ff/#fifo9" ; } 
 { { rm -f "/tmp/pash_Rx37Zq5/417af50694ad42819b8311edb72c13c3/#fifo13" ; } 
 { { rm -f "/tmp/pash_Rx37Zq5/417af50694ad42819b8311edb72c13c3/#fifo22" ; } 
 { rm -f "/tmp/pash_Rx37Zq5/417af50694ad42819b8311edb72c13c3/#fifo23" ; } ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_Rx37Zq5/9a07e1830b874c71989ee6e88728c8ff/#fifo9" ; } 
 { { mkfifo "/tmp/pash_Rx37Zq5/417af50694ad42819b8311edb72c13c3/#fifo13" ; } 
 { { mkfifo "/tmp/pash_Rx37Zq5/417af50694ad42819b8311edb72c13c3/#fifo22" ; } 
 { mkfifo "/tmp/pash_Rx37Zq5/417af50694ad42819b8311edb72c13c3/#fifo23" ; } ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ runtime/r_wrap bash -c ' tr A-Z a-z ' <"/tmp/pash_Rx37Zq5/417af50694ad42819b8311edb72c13c3/#fifo23" >"/tmp/pash_Rx37Zq5/9a07e1830b874c71989ee6e88728c8ff/#fifo9" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/r_wrap bash -c ' grep "\\(.\\).*\\1\\(.\\).*\\2\\(.\\).*\\3\\(.\\).*\\4" ' <"/tmp/pash_Rx37Zq5/9a07e1830b874c71989ee6e88728c8ff/#fifo9" >"/tmp/pash_Rx37Zq5/417af50694ad42819b8311edb72c13c3/#fifo13" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_Rx37Zq5/417af50694ad42819b8311edb72c13c3/#fifo22" -o "/tmp/pash_Rx37Zq5/417af50694ad42819b8311edb72c13c3/#fifo23" -I -m 1G -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3 aws/s3-chunk-reader-approx-correction.py "oneliners/inputs/100M.txt" "/tmp/pash_Rx37Zq5/417af50694ad42819b8311edb72c13c3/#fifo22" '[{"start": 0, "end": 3276802, "block_id": 0, "shard_id": 0}, {"start": 6553606, "end": 9830408, "block_id": 2, "shard_id": 1}, {"start": 13107212, "end": 16384014, "block_id": 4, "shard_id": 2}, {"start": 19660818, "end": 22937620, "block_id": 6, "shard_id": 3}, {"start": 26214424, "end": 29491226, "block_id": 8, "shard_id": 4}, {"start": 32768030, "end": 36044832, "block_id": 10, "shard_id": 5}, {"start": 39321636, "end": 42598438, "block_id": 12, "shard_id": 6}, {"start": 45875242, "end": 49152044, "block_id": 14, "shard_id": 7}, {"start": 52428848, "end": 55705650, "block_id": 16, "shard_id": 8}, {"start": 58982454, "end": 62259256, "block_id": 18, "shard_id": 9}, {"start": 65536060, "end": 68812862, "block_id": 20, "shard_id": 10}, {"start": 72089666, "end": 75366468, "block_id": 22, "shard_id": 11}, {"start": 78643272, "end": 81920074, "block_id": 24, "shard_id": 12}, {"start": 85196878, "end": 88473680, "block_id": 26, "shard_id": 13}, {"start": 91750484, "end": 95027286, "block_id": 28, "shard_id": 14}, {"start": 98304090, "end": 101580892, "block_id": 30, "shard_id": 15}]' shard=0 num_shards=2 job_uid=16f71f20-0d1d-4a45-9395-ec67515ee432 debug=True window_size=None chunks_per_lambda=16 write_headers=true >"/tmp/pash_Rx37Zq5/417af50694ad42819b8311edb72c13c3/#fifo22" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/pashlib-ft send*4fa45390-6e5b-4123-963e-786edd2f8586*0*1*/tmp/pash_Rx37Zq5/417af50694ad42819b8311edb72c13c3/#fifo13 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )