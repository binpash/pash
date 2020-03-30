from definitions.ir.command import *

class AltBigramGReduce(Command):
    def __init__(self, old_command, file_ids):
        assert(str(old_command.command) == "alt_bigram_aux_reduce")
        command = old_command.command
        options = self.make_options(file_ids)
        ## TODO: Generalize the model to arbitrarily many outputs to
        ## get rid of this hack, where inputs, outputs are just
        ## written as options.
        in_stream = []
        ## WARNING: This cannot be changes to be empty, because then
        ## it would not be correctly unified in the backend with the
        ## input of its downstream node.
        out_stream = [("option", 2)]
        ## TODO: When we generalize the model to have arbitrary
        ## outputs we can have them also be outputs.
        opt_indices = [("option", 0), ("option", 1)]
        category = "pure"
        super().__init__(None, command, options, in_stream, out_stream,
                         opt_indices, category)

    ## TODO: Make this atrocity prettier
    def make_options(self, file_ids):
        options = []
        for i in range(0, 3):
            fid = file_ids[i]
            if(not i == 2):
                opt = string_to_argument(fid.pipe_name())
            else:
                opt = fid
            options.append(opt)
        return options


