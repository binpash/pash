from definitions.ir.dfg_node import *

class Eager(DFGNode):
    def __init__(self, inputs, outputs, com_name, com_category, com_options = [], 
                 com_redirs = [], com_assignments=[]):
        super().__init__(inputs, outputs, com_name, com_category, 
                         com_options=com_options, 
                         com_redirs=com_redirs, 
                         com_assignments=com_assignments)

def make_eager_node(input_id, output_id, intermediate_file_id, eager_exec_path):
    com_name = Arg(string_to_argument(eager_exec_path))
    com_category = "pure"
    ## TODO: In theory the intermediate file id is also an output...
    com_options = [(2, intermediate_file_id)]
    return Eager([input_id],
                 [output_id],
                 com_name, 
                 com_category,
                 com_options=com_options)
