from definitions.ir.command import *

class SortGReduce(Command):
    def __init__(self, file_ids):
        command = string_to_argument("sort")
        input_file_ids = file_ids[:-1]
        output_file_id = file_ids[-1]
        options = [string_to_argument("-m")] + input_file_ids
        in_stream = [("option", 1), ("option", 2)]
        out_stream = ["stdout"]
        opt_indices = [("option", 0)]
        category = "pure"
        super().__init__(None, command, options, in_stream, out_stream,
                         opt_indices, category, stdout=output_file_id)

