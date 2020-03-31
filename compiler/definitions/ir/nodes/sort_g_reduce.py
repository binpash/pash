from definitions.ir.command import *

class SortGReduce(Command):
    def __init__(self, old_command, file_ids):
        assert(str(old_command.command) == "sort")
        command = old_command.command
        input_file_ids = file_ids[:-1]
        output_file_id = file_ids[-1]
        old_options = old_command.get_non_file_options()
        options = [string_to_argument("-m")] + old_options + input_file_ids
        in_stream = [("option", 1 + len(old_options)), ("option", 2 + len(old_options))]
        out_stream = ["stdout"]
        opt_indices = [("option", i) for i in range(len(old_options) + 1)]
        category = "pure"
        super().__init__(None, command, options, in_stream, out_stream,
                         opt_indices, category, stdout=output_file_id)

