
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
{ rm -f "/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo9" ; } 
 { { rm -f "/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo18" ; } 
 { rm -f "/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo19" ; } ; }
}
mkfifo_pash_fifos() {
{ mkfifo "/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo9" ; } 
 { { mkfifo "/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo18" ; } 
 { mkfifo "/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo19" ; } ; }
}
rm_pash_fifos
mkfifo_pash_fifos
pids_to_kill=""
{ cat "/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo19" | tr A-Z a-z >"/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo9" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo18" -o "/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo19" -I -m 1G -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3 aws/s3-chunk-reader-approx-correction.py "oneliners/inputs/1G.txt" "/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo18" '[{"start": 0, "end": 32767999, "block_id": 0, "shard_id": 0}, {"start": 65536000, "end": 98303999, "block_id": 2, "shard_id": 1}, {"start": 131072000, "end": 163839999, "block_id": 4, "shard_id": 2}, {"start": 196608000, "end": 229375999, "block_id": 6, "shard_id": 3}, {"start": 262144000, "end": 294911999, "block_id": 8, "shard_id": 4}, {"start": 327680000, "end": 360447999, "block_id": 10, "shard_id": 5}, {"start": 393216000, "end": 425983999, "block_id": 12, "shard_id": 6}, {"start": 458752000, "end": 491519999, "block_id": 14, "shard_id": 7}, {"start": 524288000, "end": 557055999, "block_id": 16, "shard_id": 8}, {"start": 589824000, "end": 622591999, "block_id": 18, "shard_id": 9}, {"start": 655360000, "end": 688127999, "block_id": 20, "shard_id": 10}, {"start": 720896000, "end": 753663999, "block_id": 22, "shard_id": 11}, {"start": 786432000, "end": 819199999, "block_id": 24, "shard_id": 12}, {"start": 851968000, "end": 884735999, "block_id": 26, "shard_id": 13}, {"start": 917504000, "end": 950271999, "block_id": 28, "shard_id": 14}, {"start": 983040000, "end": 1015807999, "block_id": 30, "shard_id": 15}]' shard=0 num_shards=2 job_uid=ac548bb9-6985-42d7-8eee-736062413892 debug=True window_size=None chunks_per_lambda=16 write_headers=false >"/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo18" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/pashlib-ft send*4c1a5c79-5201-459a-83ae-1a8225370842*0*1*/tmp/pash_NSonVV8/0682d5c7ec004809ad404f1434914cab/#fifo9 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )