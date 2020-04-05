from definitions.ir.command import *

class Eager(Command):
    def __init__(self, ast, command, options, in_stream, out_stream,
                 opt_indices, category, stdin=None, stdout=None):
        super().__init__(ast, command, options, in_stream, out_stream,
                         opt_indices, category, stdin=stdin, stdout=stdout)

def make_eager_node(input_file_id, output_file_id, intermediate_file_id, path):
    command = string_to_argument(path)
    options = [input_file_id, output_file_id, intermediate_file_id]
    in_stream = [("option", 0)]
    out_stream = [("option", 1)]
    ## The intermediate file is not really an option (but an output).
    ##
    ## TODO: Change that
    opt_indices = [("option", 2)]
    category = "pure"
    ## TODO: Fill the AST
    ast = None
    return Eager(ast, command, options, in_stream, out_stream,
                 opt_indices, category)

