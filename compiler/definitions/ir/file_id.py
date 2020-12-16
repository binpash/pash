from union_find import *
from ir_utils import *
from util import *

from definitions.ir.resource import *

## Note: The NULL ident is considered to be the default unknown file id
##
## TODO: WARNING: We have to make sure that a resource in our IR can
## be uniquely determined given a relative or absolute path. Actually,
## we need to make sure that expanding any variable/string in our IR,
## will always return the same result.
##
## WARNING: At the moment it is not clear what resources are saved in
## the Find(self) and in self. This might create problems down the
## road.
##
## TODO: When doing union, I have to really make both file ids point
## to the same file.
class FileId:
    def __init__(self, ident, prefix="", resource=None):
        self.ident = ident
        ## TODO: Add this as part of the resource. Ephemeral resources should
        ##       have a prefix. Normal resources have a file system id. Descriptor
        ##       resources have a number.
        self.prefix = prefix
        ## Initialize the parent
        MakeSet(self)
        self.resource=resource

    def __repr__(self):
        ## Note: Outputs the parent of the union and not the file id
        ##       itself.
        # return self.serialize()
        if (self.resource is None):
            output = "{}#file{}".format(self.prefix, Find(self).ident)
        else:
            output = "fid:{}:{}".format(Find(self).ident, self.resource)
        return output

    def serialize(self):
        # log("File id:", self.ident, Find(self).ident, self.resource, self.children)
        if (self.resource is None):
            output = "{}#file{}".format(self.prefix, Find(self).ident)
        else:
            output = "{}".format(self.resource)
        return output

    ## Serialize as an option for the JSON serialization when sent to
    ## the backend. This is needed as options can either be files or
    ## arguments, and in each case there needs to be a different
    ## serialization procedure.
    def opt_serialize(self):
        return '"{}"'.format(self.serialize())

    ## Returns a shell AST from this file identifier.
    ## TODO: Once the python libdash bindings are done we could use those instead.
    def to_ast(self):
        ## TODO: This here is supposed to identify fifos, but real fifos have a resource
        ##       but are fifos. Therefore eventually we want to have this check correctly
        ##       check if a file id refers to a pipe
        if(self.resource is None):
            string = "{}#file{}".format(self.prefix, Find(self).ident)
            ## Quote the argument
            argument = [make_kv('Q', string_to_argument(string))]
        else:
            string = "{}".format(self.resource)
            argument = string_to_argument(string)

        return argument


    ## TODO: Maybe this can be merged with serialize
    def pipe_name(self):
        assert(self.resource is None)
        output = '"#file{}"'.format(Find(self).ident)
        return output

    def set_resource(self, resource):
        ## The resource cannot be reset. A pointer can never point to
        ## more than one resource.
        assert(self.resource is None)
        self.resource = resource

    def get_resource(self):
        return self.resource

    def has_resource(self):
        return (not self.resource is None)

    def has_file_resource(self):
        return (isinstance(self.resource, FileResource))

    ## Removes a resource from an FID, making it ephemeral
    def make_ephemeral(self):
        self.resource = None

    def toFileName(self, prefix):
        output = "{}_file{}".format(prefix, Find(self).ident)
        return output

    def isNull(self):
        return self.ident == "NULL"

    ## TODO: This union-find structure is very brittle. It could break
    ## very easily and it is difficult to reason about the
    ## files. Replace this with some other structure. Maybe a map that
    ## keeps the unifications of file identifiers. Make sure that when
    ## uniting two files, their children and resources are modified
    ## accordingly.
    def union(self, other):
        Union(self, other)
        my_resource = self.get_resource()
        other_resource = Find(other).get_resource()
        ## It shouldn't be the case that both resources are not NULL
        assert(my_resource is None or
               other_resource is None or
               my_resource == other_resource)

        if (my_resource is None):
            self.set_resource(other_resource)
        elif (other_resource is None):
            Find(other).set_resource(my_resource)

    def find_fid_list(self, fids):
        parent_fids = [Find(other_fid) for other_fid in fids]
        try:
            return parent_fids.index(Find(self))
        except ValueError:
            return None

    def get_ident(self):
        return self.ident

    ##  TODO: Remove
    def flatten(self):
        if(len(self.get_children()) > 0):
            return flatten_list([child.flatten() for child in self.get_children()])
        else:
            return [self]

