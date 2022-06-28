from datatypes_new.AccessKind import AccessKind
from datatypes_new.BasicDatatypes import Flag, ArgStringType
from datatypes_new.BasicDatatypesWithIO import OptionWithIO
from datatypes_new.CommandInvocationWithIOVars import CommandInvocationWithIOVars

from annotations_utils.util_cmd_invocations import to_ast_flagoption, to_ast_operand
from definitions.ir.dfg_node import *

class DGSHTee(DFGNode):
    def __init__(self,
                 cmd_invocation_with_io_vars,
                 com_redirs=[], com_assignments=[]
                 ):
        # TODO []: default
        super().__init__(cmd_invocation_with_io_vars,
                         com_redirs=com_redirs,
                         com_assignments=com_assignments)

    # TODO: this is only needed since dgsh.sh does not comply with the XBD standard
    def to_ast(self, edges, drain_streams):
        if (drain_streams):
            raise NotImplementedError()
        else:
            redirs = self._to_ast_aux_get_redirs()
            assignments = self.com_assignments
            node = to_node_cmd_inv_with_io_vars_for_dgsh_tee(self.cmd_invocation_with_io_vars, edges, redirs, assignments)
        return node

def make_dgsh_tee_node(input_id, output_id):
    dgsh_tee_bin = os.path.join(config.PASH_TOP, config.config['runtime']['dgsh_tee_binary'])

    operand_list = [input_id, output_id]
    access_map = {output_id: AccessKind.make_stream_output(),
                  input_id: AccessKind.make_stream_input()}

    flag_option_list = [Flag("-I"),
                        Flag("-f"),
                        OptionWithIO("-b", ArgStringType(Arg(string_to_argument(str(config.config['runtime']['dgsh_buffer_size'])))))]

    cmd_inv_with_io_vars = CommandInvocationWithIOVars(
        cmd_name=dgsh_tee_bin,
        flag_option_list=flag_option_list,
        operand_list=operand_list,
        implicit_use_of_streaming_input=None,
        implicit_use_of_streaming_output=None,
        access_map=access_map)
    return DGSHTee(cmd_inv_with_io_vars)

def to_node_cmd_inv_with_io_vars_for_dgsh_tee(cmd_inv, edges, redirs, assignments):
    ast_cmd_name = string_to_argument(cmd_inv.cmd_name)
    ast_flagoptions = []
    for flagoption in cmd_inv.flag_option_list:
        ast_flagoptions += to_ast_flagoption(flagoption, edges)
    ast_operands = [to_ast_operand(operand, edges) for operand in cmd_inv.operand_list]
    # This is where it differs ... in the order
    cmd_asts = [ast_cmd_name] + ast_operands + ast_flagoptions
    # we omit stuff for stdin and stdout as we know it does not exist
    node = make_command(cmd_asts, redirections=redirs, assignments=assignments)
    return node
