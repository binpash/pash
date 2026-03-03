export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export LOCPATH=/var/task/runtime/locale
export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export RUST_BACKTRACE=1
version=$2
mkdir -p /tmp/pash_NSonVV8/ 
mkdir -p /tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/ 
rm_pash_fifos() {
{ rm -f "/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo15" ; } 
 { { rm -f "/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo20" ; } 
 { rm -f "/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo21" ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo15" ; } 
 { { mkfifo "/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo20" ; } 
 { mkfifo "/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo21" ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ sort <"/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo21" >"/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo15" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo20" -o "/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo21" -I -m 1G -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3 aws/s3-chunk-reader-approx-correction.py "oneliners/inputs/1G.txt" "/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo20" '[{"start": 32768000, "end": 65535999, "block_id": 1, "shard_id": 0}, {"start": 98304000, "end": 131071999, "block_id": 3, "shard_id": 1}, {"start": 163840000, "end": 196607999, "block_id": 5, "shard_id": 2}, {"start": 229376000, "end": 262143999, "block_id": 7, "shard_id": 3}, {"start": 294912000, "end": 327679999, "block_id": 9, "shard_id": 4}, {"start": 360448000, "end": 393215999, "block_id": 11, "shard_id": 5}, {"start": 425984000, "end": 458751999, "block_id": 13, "shard_id": 6}, {"start": 491520000, "end": 524287999, "block_id": 15, "shard_id": 7}, {"start": 557056000, "end": 589823999, "block_id": 17, "shard_id": 8}, {"start": 622592000, "end": 655359999, "block_id": 19, "shard_id": 9}, {"start": 688128000, "end": 720895999, "block_id": 21, "shard_id": 10}, {"start": 753664000, "end": 786431999, "block_id": 23, "shard_id": 11}, {"start": 819200000, "end": 851967999, "block_id": 25, "shard_id": 12}, {"start": 884736000, "end": 917503999, "block_id": 27, "shard_id": 13}, {"start": 950272000, "end": 983039999, "block_id": 29, "shard_id": 14}, {"start": 1015808000, "end": 1048576000, "block_id": 31, "shard_id": 15}]' shard=1 num_shards=2 job_uid=ac548bb9-6985-42d7-8eee-736062413892 debug=True window_size=None chunks_per_lambda=16 write_headers=false >"/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo20" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/pashlib-ft send*b04e1ae6-5eb1-43b1-95a0-7a8334174ba1*0*1*/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo15 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )
