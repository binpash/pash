hadoop jar $jarpath -files 1_map.sh,1_reduce.sh -D dfs.checksum.type=NULL -input $infile -output $output_dirs/max -mapper 1_map.sh -reducer 1_reduce.sh &
hadoop jar $jarpath -files 2_map.sh,2_reduce.sh -D dfs.checksum.type=NULL -input $infile -output $output_dirs/min -mapper 2_map.sh -reducer 2_reduce.sh &
hadoop jar $jarpath -files 3_map.sh,3_reduce.sh -D dfs.checksum.type=NULL -input $infile -output $output_dirs/average -mapper 3_map.sh -reducer 3_reduce.sh &
wait
