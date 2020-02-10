from definitions.ir.command import *

class BigramGMap(Command):
    num_outputs = 3
    def __init__(self, file_ids):
        command = string_to_argument("bigram_aux_map")
        options = [string_to_argument(fid.pipe_name()) for fid in file_ids]
        ## TODO: Generalize the model to arbitrarily many outputs to
        ## get rid of this hack, where inputs, outputs are just
        ## written as options.
        in_stream = []
        out_stream = []
        ## TODO: When we generalize the model to have arbitrary
        ## outputs we can have them also be outputs.
        opt_indices = [("option", 0), ("option", 1), ("option", 2), ("option", 3)]
        category = "stateless"
        super().__init__(None, command, options, in_stream, out_stream,
                         opt_indices, category)
