from definitions.ir.dfg_node import *
from ir_utils import *

class RemoteExec(DFGNode):
    def __init__(self, inputs, outputs, com_name, com_category,
                 com_options = [], com_redirs = [], com_assignments=[]):
        super().__init__(inputs, outputs, com_name, com_category,
                         com_options=com_options, 
                         com_redirs=com_redirs, 
                         com_assignments=com_assignments)
    
def make_remote_exec_node(node):
    remote_exec_bin = os.path.join(config.PASH_TOP, config.config['runtime']['remote_exec_binary'])
    com_name = Arg(string_to_argument(remote_exec_bin))
    ## TODO: Is it actually pure? What is it?
    com_category = "pure"
    ## At the moment we can only remote_exec a node that takes its input from stdin 
    ## and outputs to stdout. Therefore the node needs to have only one input and one output.
    inputs = node.inputs
    assert(is_single_input(inputs))

    outputs = node.outputs
    ## TODO: Would it make sense for outputs to be less than one?
    assert(len(outputs) <= 1)

    ## TODO: All arguments must be options, otherwise there must be
    ##       special handling in the wrap node2ast code. 
    # single_quote = Arg(string_to_argument("\""))
    cmd = Arg(string_to_argument(""))

    #create remote_exec argument
    # cmd.concatenate(single_quote)
    cmd.concatenate(node.com_name)
    for i, opt in node.com_options:
        cmd.concatenate(opt)
    # cmd.concatenate(single_quote)
    # command_arg = [(0, Arg([quote_arg(cmd.arg_char_list)]))]
    command_arg = [(0, cmd.wrap_in_quotes())]
    
    redirs = node.com_redirs

    assignments = node.com_assignments

    return RemoteExec(inputs,
                 outputs,
                 com_name, 
                 com_category,
                 com_options=command_arg,
                 com_redirs=redirs,
                 com_assignments=assignments)