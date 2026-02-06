from definitions.ir.dfg_node import *
from pash_annotations.datatypes.CommandInvocationWithIOVars import CommandInvocationWithIOVars, OptionWithIOVar
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
#here remote key is filename on s3 in our case 
def make_serverless_remote_pipe(local_fifo_id, is_remote_read, remote_key, output_edge=None, is_tcp=False, is_s3_lambda=False, lambda_counter=0, total_lambdas=0, byte_range=None, job_uid=None, skip_first_line=True, window_size=None, chunks_per_lambda=None, window_after_vec=None, write_headers=True):
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
        else: #TODO modify to python3.10 as python3.9 is deprecated in april 2026
            if is_s3_lambda:
                # Check for five possible modes
                import os
                use_adaptive_boundaries = os.environ.get('USE_ADAPTIVE_BOUNDARIES', 'false').lower() == 'true'
                use_dynamic_boundaries = os.environ.get('USE_DYNAMIC_BOUNDARIES', 'false').lower() == 'true'
                use_adaptive_simple = os.environ.get('USE_ADAPTIVE_SIMPLE', 'false').lower() == 'true'
                use_approx_with_correction = os.environ.get('USE_APPROX_LAMBDA_CORRECTION', 'false').lower() == 'true'
                use_smart_boundaries = os.environ.get('USE_SMART_BOUNDARIES', 'false').lower() == 'true'
                use_single_shot = os.environ.get('USE_SINGLE_SHOT', 'false').lower() == 'true'
                print("got to here")
                if use_adaptive_boundaries or use_dynamic_boundaries or use_adaptive_simple or use_approx_with_correction or use_single_shot:
                    # All three modes use same script, differentiated by window_size value
                    remote_pipe_bin = "python3 aws/s3-approx-boundaries-lambda-correction.py"
                elif use_smart_boundaries and not skip_first_line:
                    # Smart boundaries: pre-aligned, no inter-lambda communication
                    remote_pipe_bin = "python3 aws/s3-get-object-smart-and-streaming.py"
                else:
                    # Approximate boundaries: requires inter-lambda tail byte communication
                    print("Using approximate boundaries for S3 lambda reads")
                    remote_pipe_bin = "python3 aws/s3-shard-reader-streaming.py"

                # Add operands
                operand_list.append(Operand(Arg.string_to_arg(str(remote_key)))) # s3 key
                operand_list.append(local_fifo_id) #out fifo
                operand_list.append(Operand(Arg.string_to_arg(str(byte_range)))) # byte range
                operand_list.append(Operand(Arg.string_to_arg(f"shard={lambda_counter}"))) # shard
                operand_list.append(Operand(Arg.string_to_arg(f"num_shards={total_lambdas}"))) # num shards
                operand_list.append(Operand(Arg.string_to_arg(f"job_uid={job_uid}"))) # job uid
                operand_list.append(Operand(Arg.string_to_arg(f"debug=True")))

                # Add window_size operand for correction-based modes
                if use_adaptive_boundaries or use_dynamic_boundaries or use_adaptive_simple or use_approx_with_correction or use_single_shot:
                    if window_size == "adaptive":
                        # Signal adaptive mode to lambda
                        operand_list.append(Operand(Arg.string_to_arg(f"window_size=adaptive")))

                        # Pass adaptive configuration
                        target_lines = int(os.environ.get('PASH_ADAPTIVE_TARGET_LINES', '500'))
                        retry_risk = float(os.environ.get('PASH_ADAPTIVE_RETRY_RISK', '0.001'))
                        sample_kb = int(os.environ.get('PASH_ADAPTIVE_SAMPLE_KB', '256'))
                        safety_factor = float(os.environ.get('PASH_ADAPTIVE_SAFETY_FACTOR', '1.2'))

                        operand_list.append(Operand(Arg.string_to_arg(f"target_lines={target_lines}")))
                        operand_list.append(Operand(Arg.string_to_arg(f"retry_risk={retry_risk}")))
                        operand_list.append(Operand(Arg.string_to_arg(f"sample_kb={sample_kb}")))
                        operand_list.append(Operand(Arg.string_to_arg(f"safety_factor={safety_factor}")))
                    elif window_size is not None:
                        # Numeric window_size
                        operand_list.append(Operand(Arg.string_to_arg(f"window_size={window_size}")))
                    else:
                        # Dynamic mode (window_size=None)
                        operand_list.append(Operand(Arg.string_to_arg("window_size=None")))

                    if window_after_vec is not None:
                        operand_list.append(Operand(Arg.string_to_arg(f"window_after_vec={window_after_vec}")))

                # Add chunks_per_lambda operand for correction-based modes
                if (use_adaptive_boundaries or use_dynamic_boundaries or use_adaptive_simple or use_approx_with_correction or use_single_shot) and chunks_per_lambda is not None:
                    operand_list.append(Operand(Arg.string_to_arg(f"chunks_per_lambda={chunks_per_lambda}")))

                # Add write_headers operand (only if explicitly False)
                if write_headers is False:
                    operand_list.append(Operand(Arg.string_to_arg("write_headers=false")))

                # CRITICAL: Set output edge for all s3_lambda modes                                            │
                implicit_use_of_streaming_output = local_fifo_id
            else:
                remote_pipe_bin = "python3.9 aws/s3-get-object.py"
                operand_list.append(Operand(Arg.string_to_arg(str(remote_key))))
                if output_edge and (output_edge.get_resource() is not None):
                    # we need to redirect the output to some file or stdout
                    operand_list.append(Operand(Arg.string_to_arg("/dev/stdout")))
                    implicit_use_of_streaming_output = output_edge.get_ident() 
                else:
                    operand_list.append(local_fifo_id) # here we should also append byte range with the correct range
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
