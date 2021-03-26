from definitions.ir.dfg_node import *
from ir_utils import *

class RUnwrap(DFGNode):
    def __init__(self, inputs, outputs, com_name, com_category,
                 com_options = [], com_redirs = [], com_assignments=[]):
        super().__init__(inputs, outputs, com_name, com_category,
                         com_options=com_options, 
                         com_redirs=com_redirs, 
                         com_assignments=com_assignments)

def make_unwrap_node(inputs, output):
    assert(is_single_input(inputs))
    r_unwrap_bin = os.path.join(config.PASH_TOP, config.config['runtime']['r_unwrap_binary'])
    com_name = Arg(string_to_argument(r_unwrap_bin))
    com_category = "pure"
    return RUnwrap(inputs,
                   [output],
                   com_name, 
                   com_category)
