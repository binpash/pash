export PATH=$PATH:runtime
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:runtime/lib
export LOCPATH=/var/task/runtime/locale
export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export RUST_BACKTRACE=1
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
{ sort <"/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo21" >"/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo15" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/dgsh-tee -i "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo20" -o "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo21" -I -m 1G -b 5M & }
pids_to_kill="${!} ${pids_to_kill}"
{ python3 aws/s3-chunk-reader-approx-correction.py "oneliners/inputs/1M.txt" "/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo20" '[{"start": 524288, "end": 1048576, "block_id": 1, "shard_id": 0}]' shard=1 num_shards=2 job_uid=efba3ebf-ca02-4246-a4a0-504288102afc debug=True window_size=None chunks_per_lambda=1 write_headers=false >"/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo20" & }
pids_to_kill="${!} ${pids_to_kill}"
{ runtime/pashlib-ft send*53880cce-bf65-4a27-b383-c4ad58304e7c*0*1*/tmp/pash_mExGzpm/5e70c60be6cb4f048462289dd017d2ad/#fifo15 & }
pids_to_kill="${!} ${pids_to_kill}"
source runtime/wait_for_output_and_sigpipe_rest.sh ${pids_to_kill}
rm_pash_fifos
( exit "${internal_exec_status}" )