from definitions.ir.dfg_node import *

class Cat(DFGNode):
    def __init__(self, inputs, outputs, com_name, com_category, com_options = [], 
                 com_redirs = [], com_assignments=[]):
        super().__init__(inputs, outputs, com_name, com_category, 
                         com_options=com_options, 
                         com_redirs=com_redirs, 
                         com_assignments=com_assignments)

def make_cat_node(input_file_ids, output_file_id):
    raise NotImplementedError()
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
