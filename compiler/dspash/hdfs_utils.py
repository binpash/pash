from utils import get_cmd_output
from hdfs_file_data import FileData

def get_cmd_output(cmd):
    ret = subprocess.check_output(cmd, shell=True, universal_newlines=True, stderr=subprocess.PIPE)
    return ret.strip()
    
def _remove_prefix(s, prefix):
    if s.startswith(prefix):
        return s[len(prefix):]
    return s

def get_datanode_dir():
    data_dir = get_cmd_output("hdfs getconf -confKey dfs.datanode.data.dir")
    data_dir = _remove_prefix(data_dir, "file://")
    return data_dir

def get_file_data(filename):
    return FileData(filename)