from definitions.ir.arg import *
from util import *
from ir_utils import *
import socket

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
    
    def is_available_on(self, host):
        return False
    
class FileDescriptorResource(Resource):
    def __init__(self, fd):
        assert(isinstance(fd, tuple)
               and len(fd) == 2
               and fd[0] == 'fd')
        self.uri = fd

    def is_stdin(self):
        return (self.uri == ('fd', 0))

    def is_stdout(self):
        return (self.uri == ('fd', 1))

    def is_available_on(self, host):
        return True


class FileResource(Resource):
    ## The uri is the path of the file.
    def __init__(self, path):
        assert(isinstance(path, Arg))
        ## TODO: Make sure that paths are normalized
        self.uri = path
        self.data_hosts = [socket.gethostbyname(socket.gethostname())]

    def __eq__(self, other):
        if isinstance(other, FileResource):
            return self.uri == other.uri
        return False
    
    def is_available_on(self, host):
        return host in self.data_hosts


class EphemeralResource(Resource):
    def __init__(self):
        self.uri = None

    def is_available_on(self, host):
        return True
