from definitions.ir.dfg_node import *
from ir_utils import *

class RWrap(DFGNode):
    def __init__(self, inputs, outputs, com_name, com_category,
                 com_options = [], com_redirs = [], com_assignments=[], wrapped_node_name=None):
        super().__init__(inputs, outputs, com_name, com_category,
                         com_options=com_options, 
                         com_redirs=com_redirs, 
                         com_assignments=com_assignments)
        self.wrapped_node_name = wrapped_node_name
    
    ## Get the label of the node. By default, it is simply the name
    def get_dot_label(self) -> str:
        ## The name could be a full path
        name = self.com_name
        basename = os.path.basename(str(name))

        wrapped_node_name = self.wrapped_node_name
        return f'{basename}({wrapped_node_name})'

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
    single_quote = Arg(string_to_argument("\'"))
    cmd = Arg(string_to_argument(""))

    #create bash -c argument
    cmd.concatenate(single_quote)
    cmd.concatenate(node.com_name)
    for i, opt in node.com_options:
        cmd.concatenate(opt)
    cmd.concatenate(single_quote)

    wrapped_command_arg = [(1, cmd)]
    bash_command_arg = [(0, Arg(string_to_argument("bash -c")))]
    options = bash_command_arg +  wrapped_command_arg

    ## TODO: It is not clear if it is safe to just pass redirections and assignments down the line as is
    redirs = node.com_redirs
    assignments = node.com_assignments

    return RWrap(inputs,
                 outputs,
                 com_name, 
                 com_category,
                 com_options=options,
                 com_redirs=redirs,
                 com_assignments=assignments,
                 wrapped_node_name=node.com_name)
