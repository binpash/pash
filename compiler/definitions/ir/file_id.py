from union_find import *
from util import *

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
    def __init__(self, ident, resource=None, children = [], max_length = None):
        self.ident = ident
        ## Initialize the parent
        MakeSet(self)
        self.resource=resource
        self.children = children
        ## Max length shows what is the maximum possible length that
        ## this file shows to. Its use is mostly to split intermediate
        ## streams.
        self.max_length = max_length

    def __repr__(self):
        ## Note: Outputs the parent of the union and not the file id
        ##       itself.
        return self.serialize()
        # print("Repr File id:", self.ident, Find(self).ident, self.resource, self.children)
        # if (self.resource is None):
        #     if(self.max_length is None):
        #         output = "#file{}".format(Find(self).ident)
        #     else:
        #         output = "#file{}[max:{}]".format(Find(self).ident, self.max_length)
        # else:
        #     output = "#file{}({})".format(Find(self).ident, self.resource.__repr__())
        # return output

    def serialize(self):
        # print("File id:", self.ident, Find(self).ident, self.resource, self.children)
        if (self.resource is None):
            if(self.max_length is None):
                output = "#file{}".format(Find(self).ident)
            else:
                output = "#file{}[max:{}]".format(Find(self).ident, self.max_length)
        else:
            output = "{}".format(self.resource)
        return output

    ## Serialize as an option for the JSON serialization when sent to
    ## the backend. This is needed as options can either be files or
    ## arguments, and in each case there needs to be a different
    ## serialization procedure.
    def opt_serialize(self):
        return '"{}"'.format(self.serialize())

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

    ## TODO: We might need to reconstruct the parents from children,
    ## so we might have to add a parent field in file ids.
    def set_children(self, children):
        assert(self.children == [])
        self.children = children

    ## TODO: DO NOT FORGET: These children functions should be removed
    ## after the model changes. There are no children, just many
    ## inputs or outputs, and how each command interprets them depends
    ## on it.
    def unset_children(self):
        self.children = []

    def get_children(self):
        return self.children

    def split_resource(self, times, batch_size, fileIdGen):
        if(not self.resource is None):
            ## This works as expected if the file points to a resource
            resources = self.resource.split_resource(times, batch_size)
            split_file_ids = [create_file_id_for_resource(resource, fileIdGen)
                              for resource in resources]
            self.set_children(split_file_ids)
        else:
            ## If the file doesn't point to a resource (meaning that
            ## it is an intermediate file), then we just restrict its
            ## max length, and the implementation should know to only
            ## pass so many lines in this file.
            split_file_ids = [create_split_file_id(batch_size, fileIdGen)
                              for i in range(times)]
            split_file_ids[-1].set_max_length(None)
            self.set_children(split_file_ids)

    ## This must be used by the implementation to only transfer
    ## max_length lines in this file.
    def set_max_length(self, max_length):
        self.max_length = max_length

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

    def flatten(self):
        if(len(self.get_children()) > 0):
            return flatten_list([child.flatten() for child in self.get_children()])
        else:
            return [self]

