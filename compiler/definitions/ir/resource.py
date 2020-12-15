from util import *

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
