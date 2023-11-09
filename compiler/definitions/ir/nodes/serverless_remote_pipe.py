from definitions.ir.dfg_node import *
from pash_annotations.datatypes.CommandInvocationWithIOVars import CommandInvocationWithIOVars
from pash_annotations.datatypes.AccessKind import make_stream_input, make_stream_output
from pash_annotations.datatypes.BasicDatatypes import Operand

class ServerlessRemotePipe(DFGNode):
    def __init__(self,
                 cmd_invocation_with_io_vars,
                 com_redirs=[],
                 com_assignments=[],
                 parallelizer_list=None,
                 cmd_related_properties=None):
        super().__init__(cmd_invocation_with_io_vars=cmd_invocation_with_io_vars,
                         com_redirs=com_redirs,
                         com_assignments=com_assignments,
                         parallelizer_list=parallelizer_list,
                         cmd_related_properties=cmd_related_properties)


def make_serverless_remote_pipe(input_id, output_id, is_remote_read, key, out_resource=None, last_subgraph=False):
    operand_list = [Operand(Arg.string_to_arg(str(key)))]
    access_map = {}
    implicit_use_of_streaming_output = None
    if input_id is not None:
        access_map[input_id] = make_stream_input()
    if output_id is not None:
        access_map[output_id] = make_stream_output()
    if is_remote_read:
        remote_pipe_bin = os.path.join('aws','remote-read.py')
        if isinstance(out_resource, FileDescriptorResource):
            operand_list.append(Operand(Arg.string_to_arg("/dev/stdout")))
            implicit_use_of_streaming_output = output_id # avoid node not found err
        else:
            operand_list.append(output_id)
    else:
        remote_pipe_bin = os.path.join('aws','remote-write.py')
        operand_list.append(input_id)
        implicit_use_of_streaming_output = output_id
        if last_subgraph:
            operand_list.append(Operand(Arg.string_to_arg("1")))
        else:
            operand_list.append(Operand(Arg.string_to_arg("0")))
    
    cmd_inv_with_io_vars = CommandInvocationWithIOVars(
        cmd_name="python "+remote_pipe_bin,
        flag_option_list=[],
        operand_list=operand_list,
        implicit_use_of_streaming_input=None,
        implicit_use_of_streaming_output=implicit_use_of_streaming_output,
        access_map=access_map)
    return ServerlessRemotePipe(cmd_inv_with_io_vars)
