from definitions.ir.command import *

class Cat(Command):
    def __init__(self, ast, command, options, in_stream, out_stream,
                 opt_indices, category, stdin=None, stdout=None, redirections=[]):
        super().__init__(ast, command, options, in_stream, out_stream,
                         opt_indices, category, stdin=stdin, stdout=stdout, redirections=redirections)

def make_cat_node(input_file_ids, output_file_id):
    command = string_to_argument("cat")
    options = input_file_ids
    in_stream = [("option", i)  for i in range(len(input_file_ids))]
    out_stream = ["stdout"]
    stdout = output_file_id
    opt_indices = []
    category = "stateless"
    ## TODO: Fill the AST
    ast = None
    return Cat(ast, command, options, in_stream, out_stream,
               opt_indices, category, stdout=stdout)

