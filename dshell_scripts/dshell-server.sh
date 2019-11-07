#!/bin/bash

PORT=$1
if [ -z $PORT ]
then
	echo 'Error! => Specify the server port'
	exit 1
fi

# HDFS check
HDFS_PATH=~/hadoop-3.2.0/sbin/start-dfs.sh
if [ ! -f "$HDFS_PATH" ] 
then
	echo 'Error! => HDFS not found on current machine'
	exit 1
fi

# kill the process that uses the $PORT if it's used
process_id=$(lsof -i:$PORT -t)
if [ ! -z "$process_id" ] 
then
    kill -9 $process_id
fi

echo "Starting dshell server..."
sudo java -cp dshell.jar dshell.core.worker.DistributedTask 4000
