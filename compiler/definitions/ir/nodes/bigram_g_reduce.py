from definitions.ir.command import *

class BigramGReduce(Command):
    def __init__(self, old_command, file_ids):
        assert(str(old_command.command) == "bigrams_aux")
        command = string_to_argument("bigram_aux_reduce")
        options = self.make_options(file_ids)
        ## TODO: Generalize the model to arbitrarily many outputs to
        ## get rid of this hack, where inputs, outputs are just
        ## written as options.
        in_stream = [("option", i) for i in range(6)]
        ## WARNING: This cannot be changes to be empty, because then
        ## it would not be correctly unified in the backend with the
        ## input of its downstream node.
        out_stream = [("option", 6)]
        ## TODO: When we generalize the model to have arbitrary
        ## outputs we can have them also be outputs.
        opt_indices = [("option", i) for i in range(7,9)]
        category = "pure"
        super().__init__(None, command, options, in_stream, out_stream,
                         opt_indices, category)

    ## TODO: Make this atrocity prettier
    def make_options(self, file_ids):
        options = []
        for i in range(0, 9):
            fid = file_ids[i]
            if(i > 6):
                opt = string_to_argument(fid.pipe_name())
            else:
                opt = fid
            options.append(opt)
        return options


