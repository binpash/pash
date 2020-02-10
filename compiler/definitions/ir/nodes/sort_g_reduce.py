from definitions.ir.command import *

class SortGReduce(Command):
    def __init__(self, file_ids):
        command = string_to_argument("sort")
        input_file_ids = file_ids[:-1]
        output_file_id = file_ids[-1]
        input_file_id_opts = [string_to_argument(fid.pipe_name()) for fid in input_file_ids]
        options = [string_to_argument("-m")] + input_file_id_opts
        ## Question: Can putting an empty input stream be a problem?
        ## If we do that, then this node will be a source in the
        ## dataflow graph. Ideally we want to generalize to arbitrary
        ## inputs. This also applies to the one above.
        in_stream = []
        out_stream = ["stdout"]
        ## At the moment none of the $IN1, $IN2 are considered inputs
        ##
        ## TODO: When we generalize the model to have arbitrary
        ## inputs we can have them also be inputs.
        opt_indices = [("option", 0), ("option", 1), ("option", 2)]
        category = "pure"
        super().__init__(None, command, options, in_stream, out_stream,
                         opt_indices, category, stdout=output_file_id)

