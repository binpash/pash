## TODO: Resources should probably be more elaborate than just a
## string and a line range. They could be URLs, and possibly other things.
class Resource:
    def __init__(self, uri, range = [0, "inf"]):
        self.uri = uri
        self.range = range

    def __repr__(self):
        if self.range[0] == 0 and self.range[1] == "inf":
            output = str(self.uri)
        else:
            output = "{}[{}:{}]".format(self.uri, self.range[0], self.range[1])
        return output

    ## This function splits a "splittable" resource (at the moment
    ## only files are resources so they are all splittable)
    def split_resource(self, times, batch_size):
        assert(times > 0 and batch_size > 0)
        init = self.range[0]
        resources = [Resource(self.uri, [i * batch_size + init, (i+1) * batch_size + init])
                     for i in range(times)]
        ## Set the end to be equal to the old end
        resources[-1].range[1] = self.range[1]
        return resources

    ## This returns the length of the resource or None if it is
    ## infinite.
    def get_length(self):
        if(self.range[1] == "inf"):
            return None
        else:
            return (self.range[1] - self.range[0])
