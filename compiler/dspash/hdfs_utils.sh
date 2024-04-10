# Helper functions and env for hdfs support
# source $PASH_TOP/compiler/dspash/hdfs_utils.sh

datanode_dir=$(hdfs getconf -confKey dfs.datanode.data.dir) 
export HDFS_DATANODE_DIR=${datanode_dir#"file://"} # removes file:// prefix

function get_hdfs_block_path() {
    dnodeName="$1"
    blockID="$2"
    find "$HDFS_DATANODE_DIR/current/$dnodeName" -name "$blockID"
}

export -f get_hdfs_block_path