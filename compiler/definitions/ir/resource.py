from definitions.ir.arg import *
from util import *
from ir_utils import *

## TODO: Resources should probably be more elaborate than just a
## string and a line range. They could be URLs, and possibly other things.

## TODO: Think if we can have any optimizations if we know the size of a resource.
class Resource:
    def __init__(self, uri):
        self.uri = uri

    def __repr__(self):
        output = str(self.uri)
        return output

    ## TODO: Make this check properly
    def __eq__(self, other):
        log("          Self URI:", self, type(self.uri))
        log("          Other URI:", other, type(other))
        if isinstance(other, Resource):
            log("          Other URI:", other, type(other.uri))
            return self.uri == other.uri
        return False
    
    def is_stdin(self):
        return (self.uri == ('fd', 0))

    def is_stdout(self):
        return (self.uri == ('fd', 1))

class FileDescriptorResource(Resource):
    def __init__(self, fd):
        assert(isinstance(fd, tuple)
               and len(fd) == 2
               and fd[0] == 'fd')
        self.uri = fd

class FileResource(Resource):
    ## The uri is the path of the file.
    def __init__(self, path):
        assert(isinstance(path, Arg))
        ## TODO: Make sure that paths are normalized
        self.uri = path

class EphemeralResource(Resource):
    def __init__(self):
        self.uri = None
