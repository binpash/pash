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

def make_serverless_remote_pipe(local_fifo_id, is_remote_read, remote_key, output_edge=None, is_tcp=False):
    """
    Generate a dfg node for serverless remote communication, now we handle these cases
    1. Remote read from tcp communication
    2. Remote read from s3
        (1) Read from s3 to fifo
        (2) Read from s3 to stdout or other fd
    3. Remote write to tcp communication
    4. Remote write to s3

    Note: for 1/2, output_edge will be stdout (just to make the sink node has an output edge manually)
    """
    operand_list = []
    access_map = {}
    implicit_use_of_streaming_input = implicit_use_of_streaming_output = None
    if local_fifo_id:
        access_map[local_fifo_id] = make_stream_input()
    if output_edge:
        access_map[output_edge.get_ident()] = make_stream_output()
    if is_remote_read:
        if is_tcp:
            remote_pipe_bin = "/opt/pashlib"
            operand_list.append(Operand(Arg.string_to_arg("recv "+str(remote_key)+" 1 0")))
            implicit_use_of_streaming_output = local_fifo_id
        else:
            remote_pipe_bin = "python3.9 aws/s3-get-object.py"
            operand_list.append(Operand(Arg.string_to_arg(str(remote_key))))
            if output_edge and (output_edge.get_resource() is not None):
                # we need to redirect the output to some file or stdout
                operand_list.append(Operand(Arg.string_to_arg("/dev/stdout")))
                implicit_use_of_streaming_output = output_edge.get_ident() 
            else:
                operand_list.append(local_fifo_id)
    else:
        implicit_use_of_streaming_output = output_edge.get_ident() # avoid node not found err
        if is_tcp:
            remote_pipe_bin = "/opt/pashlib"
            operand_list.append(Operand(Arg.string_to_arg("send "+str(remote_key)+" 0 1")))
            implicit_use_of_streaming_input = local_fifo_id
        else:
            remote_pipe_bin = "python3.9 aws/s3-put-object.py"
            operand_list.append(Operand(Arg.string_to_arg(str(remote_key))))
            operand_list.append(local_fifo_id)
            operand_list.append(Arg.string_to_arg("$1"))
    cmd_inv_with_io_vars = CommandInvocationWithIOVars(
        cmd_name=remote_pipe_bin,
        flag_option_list=[],
        operand_list=operand_list,
        implicit_use_of_streaming_input=implicit_use_of_streaming_input,
        implicit_use_of_streaming_output=implicit_use_of_streaming_output,
        access_map=access_map)
    return ServerlessRemotePipe(cmd_inv_with_io_vars)


def make_serverless_remote_pipe_one_proc(list_of_arg):
    operand_list = []

    for arg in list_of_arg:
        operand_list.append(Operand(Arg.string_to_arg(arg)))

    cmd_inv_with_io_vars = CommandInvocationWithIOVars(
        cmd_name="/opt/pashlib",
        flag_option_list=[],
        operand_list=operand_list,
        access_map={},
        implicit_use_of_streaming_input=None,
        implicit_use_of_streaming_output=None
        )
    

    return ServerlessRemotePipe(cmd_inv_with_io_vars)

def make_serverless_ingate(data_type, list_of_arg):
    if data_type == "batch":
        runtime = "runtime/ingate_stateless"
    else:
        runtime = "runtime/ingate_stateful"
    
    operand_list = []

    for arg in list_of_arg:
        operand_list.append(Operand(Arg.string_to_arg(arg)))

    cmd_inv_with_io_vars = CommandInvocationWithIOVars(
        cmd_name=runtime,
        flag_option_list=[],
        operand_list=operand_list,
        access_map={},
        implicit_use_of_streaming_input=None,
        implicit_use_of_streaming_output=None
        )

    return ServerlessRemotePipe(cmd_inv_with_io_vars)

def make_serverless_outgate(data_type, list_of_arg):
    if data_type == "batch":
        runtime = "runtime/outgate_stateless"
    else:
        runtime = "runtime/outgate_stateful"

    operand_list = []

    for arg in list_of_arg:
        operand_list.append(Operand(Arg.string_to_arg(arg)))

    cmd_inv_with_io_vars = CommandInvocationWithIOVars(
        cmd_name=runtime,
        flag_option_list=[],
        operand_list=operand_list,
        access_map={},
        implicit_use_of_streaming_input=None,
        implicit_use_of_streaming_output=None
        )

    return ServerlessRemotePipe(cmd_inv_with_io_vars)