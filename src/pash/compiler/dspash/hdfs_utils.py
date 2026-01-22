from dspash.hdfs_file_data import get_hdfs_file_data, FileData, HDFSFileConfig
from typing import List, Tuple


def get_cmd_output(cmd: str):
    ret = subprocess.check_output(
        cmd, shell=True, universal_newlines=True, stderr=subprocess.PIPE
    )
    return ret.strip()


def _remove_prefix(s: str, prefix: str) -> str:
    if s.startswith(prefix):
        return s[len(prefix) :]
    return s


def get_datanode_dir() -> str:
    data_dir = get_cmd_output("hdfs getconf -confKey dfs.datanode.data.dir")
    data_dir = _remove_prefix(data_dir, "file://")
    return data_dir


def get_file_data(filename: str) -> FileData:
    return get_hdfs_file_data(filename)


def get_file_config(filename: str) -> HDFSFileConfig:
    filedata = get_file_data(filename)
    return HDFSFileConfig(filedata)
