from definitions.ir.dfg_node import *
from ir_utils import *

class RWrap(DFGNode):
    def __init__(self, inputs, outputs, com_name, com_category,
                 com_options = [], com_redirs = [], com_assignments=[]):
        super().__init__(inputs, outputs, com_name, com_category,
                         com_options=com_options, 
                         com_redirs=com_redirs, 
                         com_assignments=com_assignments)

def wrap_node(node):
    r_wrap_bin = os.path.join(config.PASH_TOP, config.config['runtime']['r_wrap_binary'])
    com_name = Arg(string_to_argument(r_wrap_bin))
    ## TODO: Is it actually pure? What is it?
    com_category = "pure"
    ## At the moment we can only wrap a node that takes its input from stdin 
    ## and outputs to stdout. Therefore the node needs to have only one input and one output.
    inputs = node.inputs
    assert(is_single_input(inputs))

    outputs = node.outputs
    ## TODO: Would it make sense for outputs to be less than one?
    assert(len(outputs) <= 1)

    ## TODO: For now we can only wrap stateless commands
    assert(node.com_category == "stateless")

    ## TODO: All arguments must be options, otherwise there must be
    ##       special handling in the wrap node2ast code. 
    old_options_transposed = [(i+1, opt) for i, opt in node.com_options]
    options = [(0, node.com_name)] + old_options_transposed

    ## TODO: It is not clear if it is safe to just pass redirections and assignments down the line as is
    redirs = node.com_redirs
    assignments = node.com_assignments

    return RWrap(inputs,
                 outputs,
                 com_name, 
                 com_category,
                 com_options=options,
                 com_redirs=redirs,
                 com_assignments=assignments)
