import json
from collections import namedtuple
from threading import Event, Thread
from typing import Callable, List, Set

import requests

# if you are running this outside of the docker container
# you may want to change this to localhost for testing
host = "namenode"
port = "9870"

daemon_quit = Event()


HDFSBlock = namedtuple("HDFSBlock", ["path", "hosts"])


# naming of this class and it's functionality is not ideal ¯\_(ツ)_/¯
# however, this class has hard to miss dependencies so it's hard to modify
# for example I was thinking about removing the dumps() method as I was thinking
# this class is only written but not read. However, it seems there may be go client
# code that reads it. See dfs_split_reader.go
class HDFSFileConfig:
    def __init__(self, blocks: List[List[str]]):
        self.blocks: List[HDFSBlock] = []
        for inner in blocks:
            # get_hdfs_block_path is a helper function defined in hdfs_utils.sh
            # it takes two arguments: directory name and block id and returns the path of the block
            # however here, path is not an exact path but a command that will be invoked on workers
            path = f"$(get_hdfs_block_path {inner[0]} {inner[1]})"
            hosts = inner[2:]
            self.blocks.append(HDFSBlock(path, hosts))

    def _serialize(self):
        data = {"blocks": []}
        for path, hosts in self.blocks:
            data["blocks"].append({"path": path, "hosts": hosts})
        return data

    def dumps(self):
        data = self._serialize()
        return json.dumps(data)

    def dump(self, filepath):
        data = self._serialize()
        with open(filepath, "w") as f:
            json.dump(data, f)

    def __eq__(self, __o: object) -> bool:
        if not isinstance(__o, HDFSFileConfig):
            return False
        return self.blocks == __o.blocks


def file_to_blocks(filepath: str) -> List[List[str]]:
    """
    Takes an hdfs file path as an input and returns a list of inner lists.
    For each inner list, following are true:
        - corresponds to a block
        - first element is the directory name used by hdfs_utils.sh
        - second element is the block id
        - rest of the elements are the ip addresses of the datanodes that have the block

    Example output:
    [['BP-68286741-172.20.0.2-1700503545710', 'blk_1073741830', '172.22.0.3', '172.22.0.4', '172.22.0.7'],
     ['BP-68286741-172.20.0.2-1700503545710', 'blk_1073741831', '172.22.0.3', '172.22.0.5', '172.22.0.7'],
     ['BP-68286741-172.20.0.2-1700503545710', 'blk_1073741832', '172.22.0.3', '172.22.0.5', '172.22.0.4'],
     ['BP-68286741-172.20.0.2-1700503545710', 'blk_1073741833', '172.22.0.5', '172.22.0.6', '172.22.0.7']]
    """
    outer = []

    url = f"http://{host}:{port}/fsck?ugi=root&files=1&blocks=1&locations=1&path={filepath}"
    r = requests.get(url=url)

    save_blocks = False
    for line in r.text.splitlines():
        if line.startswith(filepath):
            size = int(line.split()[1])
            assert size > 0
            save_blocks = True
            continue

        if save_blocks:
            if len(line) == 0:
                break

            space_ix = line.find(" ")
            semi_ix = line.find(":")
            under_ix = line.find("_", semi_ix + 5)

            dir_name = line[space_ix + 1 : semi_ix]
            block_id = line[semi_ix + 1 : under_ix]

            inner = []
            inner.append(dir_name)
            inner.append(block_id)

            after = 0
            while True:
                # len("DatanodeInfoWithStorage") + 1 = 24
                ip_ix = line.find("DatanodeInfoWithStorage", after) + 24

                # -1 + 24 = 23
                if ip_ix == 23:
                    break

                comma_ix = line.find(",", ip_ix)
                ip_addr = line[ip_ix:comma_ix]
                hostname = ip_addr.split(':')[0]
                after = comma_ix

                inner.append(hostname)

            outer.append(inner)

    return outer


def block_to_nodes(block_id: str) -> List[str]:
    """
    Takes a block id as an input and returns a list.
    First element of the list is the hdfs file path this block belongs to.
    Rest of the elements are the ip addresses of the datanodes that have the block.

    Example:
    input: blk_1073741830
    output: ['/500mib-file', '172.22.0.3', '172.22.0.4', '172.22.0.7']
    """
    res = []

    url = f"http://{host}:{port}/fsck?ugi=root&blockId={block_id}+&path=%2F"
    t = requests.get(url=url).text

    # len("Block belongs to: ") = 18
    file_ix_start = t.find("Block belongs to: ") + 18
    file_ix_end = t.find("\n", file_ix_start)

    filepath = t[file_ix_start:file_ix_end]
    res.append(filepath)

    all_blocks = file_to_blocks(filepath)
    for block in all_blocks:
        if block[1] == block_id:
            for addr in block[2:]:
                res.append(addr.split(':')[0])
            break

    return res


def get_live_nodes():
    """
    Returns a dictionary where keys are the ip addresses of the datanodes and values are some related information.
    Please be careful as the keys can contain hostnames.

    Example output:
    {
        "c107c1d2c0f0:9866": {
            "infoAddr": "172.22.0.5:9864",
            "infoSecureAddr": "172.22.0.5:0",
            "xferaddr": "172.22.0.5:9866",
            "lastContact": 0,
            "usedSpace": 393220096,
            "adminState": "In Service",
            "nonDfsUsedSpace": 16368644096,
            "capacity": 1081101176832,
            "numBlocks": 8,
            "version": "3.2.2",
            "used": 393220096,
            "remaining": 1009346957312,
            "blockScheduled": 0,
            "blockPoolUsed": 393220096,
            "blockPoolUsedPercent": 0.03637218,
            "volfails": 0,
            "lastBlockReport": 136
        },
        "15d32bc24bfd:9866": {
            "infoAddr": "172.22.0.3:9864",
            "infoSecureAddr": "172.22.0.3:0",
            "xferaddr": "172.22.0.3:9866",
            ...
        },
        "16489dccb5b2:9866": {
            ...
        },
        "27c75d6187d8:9866": {
            ...
        },
        "5783c1a1a370:9866": {
            ...
        }
    }
    """
    query = "Hadoop:service=NameNode,name=NameNodeInfo"
    url = f"http://{host}:{port}/jmx?qry={query}"
    r = requests.get(url)

    return json.loads(json.loads(r.text)["beans"][0]["LiveNodes"])


def get_active_node_addresses(lastContactThreshold: int):
    """
    Returns a list of ip addresses of the datanodes that are active.
    If a datanode has not contacted the namenode in the last lastContactThreshold seconds, it is considered inactive.

    Example output:
    ['172.25.0.5', '172.25.0.4', '172.25.0.3']
    """
    live_nodes = get_live_nodes()
    return {
        v["infoAddr"].split(":")[0]
        for v in live_nodes.values()
        if v["lastContact"] < lastContactThreshold
    }


def __hdfs_daemon(
    lastContactThreshold: int,
    initial_addresses: Set[str],
    func_added: Callable[[str], None],
    func_removed: Callable[[str], None],
):
    daemon_state = initial_addresses
    while not daemon_quit.is_set():
        daemon_quit.wait(1)
        new_deamon_state = get_active_node_addresses(lastContactThreshold)
        if new_deamon_state != daemon_state:
            added_addresses = new_deamon_state - daemon_state
            removed_addresses = daemon_state - new_deamon_state

            for addr in added_addresses:
                if func_added is not None:
                    func_added(addr)

            for addr in removed_addresses:
                if func_removed is not None:
                    func_removed(addr)

        daemon_state = new_deamon_state


def start_hdfs_daemon(
    lastContactThreshold: int,
    initial_addresses: Set[str],
    func_added: Callable[[str], None],
    func_removed: Callable[[str], None],
):
    """
    Starts a daemon that checks the active datanodes every second.
    If a datanode becomes active or inactive, it calls the corresponding function.
    Active datanodes are the ones that have contacted the namenode in the last lastContactThreshold seconds.
    initial_addresses is a set of ip addresses without ports of the datanodes that are active at the beginning.
    """
    Thread(
        target=__hdfs_daemon,
        args=(lastContactThreshold, initial_addresses, func_added, func_removed),
    ).start()


def stop_hdfs_daemon():
    daemon_quit.set()


def get_file_config(filepath: str) -> HDFSFileConfig:
    # Workaround included quotation marks when cat is called with this notation"${IN}"
    # TODO: this should be fixed somewhere higher in the stack
    filepath = filepath.lstrip('"').rstrip('"')
    blocks = file_to_blocks(filepath)
    return HDFSFileConfig(blocks)


# used for testing
if __name__ == "__main__":
    # print(file_to_blocks("/README.md"))
    # print(json.dumps(get_live_nodes(), indent=4))
    # print(get_active_node_addresses())
    start_hdfs_daemon(10, set(), None, None)
    # print(file_to_blocks("/500mib-file"))
    # print(block_to_nodes("blk_1073741830"))
