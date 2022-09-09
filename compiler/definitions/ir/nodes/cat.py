from datatypes_new.CommandInvocationWithIOVars import CommandInvocationWithIOVars
from definitions.ir.dfg_node import DFGNode

class Cat(DFGNode):
    def __init__(self, inputs, outputs, com_name, com_category,
                 com_options = [], com_redirs = [], com_assignments=[],
                 # BEGIN ANNO
                 flag_option_list = None,
                 positional_config_list = None,
                 positional_input_list = None,
                 positional_output_list = None,
                 implicit_use_of_stdin = None,
                 implicit_use_of_stdout = None,
                 parallelizer_list = None,
                 cmd_related_properties = None
                 # END ANNO
                 ):
        assert(False)
        assert(str(com_name) == "cat")
        assert(com_category == "stateless")
        super().__init__(inputs, outputs, com_name, com_category,
                         com_options=com_options,
                         flag_option_list=flag_option_list,
                         com_redirs=com_redirs,
                         com_assignments=com_assignments,
                         positional_config_list=positional_config_list,
                         positional_input_list=positional_input_list,
                         positional_output_list=positional_output_list,
                         implicit_use_of_stdin=implicit_use_of_stdin,
                         implicit_use_of_stdout=implicit_use_of_stdout,
                         parallelizer_list=parallelizer_list,
                         cmd_related_properties=cmd_related_properties
                         )

def make_cat_node(inputs, output):
    cmd_inv_cat = CommandInvocationWithIOVars.make_cat_command_invocation_with_io_vars(inputs, output)
    return DFGNode.make_simple_dfg_node_from_cmd_inv_with_io_vars(cmd_inv_cat)
