import os
import subprocess
import sys


class FileData(object):
    def __init__(self):
        self.blocknames = []
        self.dnodenames = []
        self.machines = []
        self.size = 0
        # self.dnodepath = subprocess.check_output("hdfs getconf -confKey dfs.datanode.data.dir", shell=True).decode("utf-8").strip("\n")
        # if (self.dnodepath.startswith("file://")):
        #    self.dnodepath = self.dnodepath[len("file://"):]

    def paths(self):
        assert len(self.blocknames) != 0
        assert len(self.dnodenames) != 0
        assert self.size > 0
        filepaths = []
        for i in range(len(self.blocknames)):
            filepaths.append(
                os.path.join(
                    "current",
                    self.dnodenames[i],
                    "current",
                    "finalized",
                    "subdir0",
                    "subdir0",
                    self.blocknames[i],
                )
            )
        return filepaths


def getIPs(raw):
    rawparts = raw.split(" ")
    ips = []
    for part in rawparts:
        index = part.find("DatanodeInfoWithStorage")
        ips.append(part[index + len("DatanodeInfoWithStorage") + 1 : part.find(",")])
    return ips


def getinfo(filename):
    info = FileData()
    log = subprocess.check_output(
        "hdfs fsck {0} -files -blocks -locations".format(filename), shell=True
    )
    count = 0
    for line in log.splitlines():
        wordarr = line.split()
        if len(wordarr) > 0 and wordarr[0].decode("utf-8") == filename and count == 0:
            info.size = int(wordarr[1])
            count += 1
        elif (
            len(wordarr) > 0
            and count > 0
            and wordarr[0][:-1].isdigit()
            and int(wordarr[0][:-1]) == count - 1
        ):
            count += 1
            rawinfo = wordarr[1].decode("utf-8").split(":")
            info.blocknames.append(rawinfo[1][0 : rawinfo[1].rfind("_")])
            info.dnodenames.append(rawinfo[0])
            stline = line.decode("utf-8")
            info.machines.append(
                getIPs(stline[stline.find("DatanodeInfoWithStorage") - 1 :])
            )
    assert len(info.blocknames) != 0
    assert len(info.dnodenames) != 0
    assert info.size > 0
    return info


if __name__ == "__main__":
    assert len(sys.argv) == 2
    filename = sys.argv[1]
    info = getinfo(filename)
    print("Size = ", info.size)
    paths = info.paths()
    for i in range(len(paths)):
        print("Machines = ", info.machines[i])
        print("Block = ", paths[i])
