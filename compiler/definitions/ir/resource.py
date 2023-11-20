import socket

from definitions.ir.arg import *
from util import *
from shell_ast.ast_util import *

## TODO: Resources should probably be more elaborate than just a
## string and a line range. They could be URLs, and possibly other things.


## TODO: Think if we can have any optimizations if we know the size of a resource.
class Resource:
    def __init__(self, uri):
        self.uri = uri

    def __repr__(self):
        output = str(self.uri)
        return output

    def is_stdin(self):
        return False

    def is_stdout(self):
        return False

    ## TODO: Make this check properly
    def __eq__(self, other):
        if isinstance(other, Resource):
            return self.uri == other.uri
        return False


class FileDescriptorResource(Resource):
    def __init__(self, fd):
        assert isinstance(fd, tuple) and len(fd) == 2 and fd[0] == "fd"
        self.uri = fd

    def is_stdin(self):
        return self.uri == ("fd", 0)

    def is_stdout(self):
        return self.uri == ("fd", 1)


class FileResource(Resource):
    ## The uri is the path of the file.
    def __init__(self, path):
        log("class of path", type(path))
        assert isinstance(path, Arg)
        ## TODO: Make sure that paths are normalized
        self.uri = path

    def __eq__(self, other):
        if isinstance(other, FileResource):
            return self.uri == other.uri
        return False


class TemporaryFileResource(Resource):
    def __init__(self):
        self.uri = None


# A FIFO.
class EphemeralResource(Resource):
    def __init__(self):
        self.uri = None


class RemoteFileResource(Resource):
    def __init__(self):
        raise NotImplementedError("RemoteFileResource is an interface")

    def __str__(self):
        raise NotImplementedError("RemoteFileResource is an interface")

    def is_available_on(self, host):
        raise NotImplementedError("RemoteFileResource is an interface")

    def _normalize_addr(self, addr):
        """
        Helper, removes port number if supplied in address and normalizes host address
        Example1: datanode1 -> 172.18.0.3
        Example2: 172.18.0.3:5555 -> 172.18.0.3
        """
        host = addr.rsplit(":", 1)[0]
        normalized_host = socket.gethostbyaddr(host)[2][0]
        return normalized_host


class HDFSFileResource(RemoteFileResource):
    ## The uri is the path of the file.
    def __init__(self, uri, resource_hosts):
        """
        Params:
            uri: Usually the path to the file. The path doesn't include the top directory
            which is different between hosts. The str function adds the prefix $HDFS_DATANODE_DIR/
            which should be defined on host machine worker environment.
            resource_hosts: the addresses of all the machines containing
            the resource.
        """
        self.uri = uri
        self.hosts = list(map(self._normalize_addr, resource_hosts))

    def __eq__(self, other):
        if isinstance(other, HDFSFileResource):
            return self.uri == other.uri
        return False

    def is_available_on(self, host):
        return host in self.hosts

    def __repr__(self):
        return f"hdfs://{self.uri}"

    def __str__(self):
        return f"$HDFS_DATANODE_DIR/{self.uri}"


# DFS logical split resource
class DFSSplitResource(RemoteFileResource):
    def __init__(self, config, config_path, split_num, hosts):
        self.config = config
        self.config_path = config_path
        self.split_num = split_num
        self.hosts = list(map(self._normalize_addr, hosts))

    def is_available_on(self, host):
        return host in self.hosts

    def set_config_path(self, config_path):
        self.config_path = config_path

    def __str__(self):
        return self.config_path
