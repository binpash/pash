from definitions.ir.command import *

class BigramGMap(Command):
    num_outputs = config.bigram_g_map_num_outputs
    def __init__(self, file_ids):
        assert(self.num_outputs == 3)
        command = string_to_argument("bigram_aux_map")
        options = self.make_options(file_ids)
        # options = [string_to_argument(fid.pipe_name()) for fid in file_ids]
        ## TODO: Generalize the model to arbitrarily many outputs to
        ## get rid of this hack, where inputs, outputs are just
        ## written as options.
        in_stream = [("option", 0)]
        out_stream = []
        ## TODO: When we generalize the model to have arbitrary
        ## outputs we can have them also be outputs.
        opt_indices = [("option", 1), ("option", 2), ("option", 3)]
        category = "stateless"
        super().__init__(None, command, options, in_stream, out_stream,
                         opt_indices, category)

    ## TODO: Make this atrocity prettier
    def make_options(self, file_ids):
        options = []
        for i in range(4):
            fid = file_ids[i]
            if(not i == 0):
                opt = string_to_argument(fid.pipe_name())
            else:
                opt = fid
            options.append(opt)
        return options
