from pash_annotations.datatypes.AccessKind import make_stream_output, make_stream_input
from pash_annotations.datatypes.BasicDatatypes import Flag, ArgStringType
from pash_annotations.datatypes.BasicDatatypesWithIO import OptionWithIO
from pash_annotations.datatypes.CommandInvocationWithIOVars import CommandInvocationWithIOVars

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

def make_dgsh_tee_node(input_id, output_id, is_s3=False):
    dgsh_tee_bin = os.path.join(config.PASH_TOP, config.config['runtime']['dgsh_tee_binary'])
    if config.pash_args.serverless_exec:
        dgsh_tee_bin = config.config['runtime']['dgsh_tee_binary']

    access_map = {output_id: make_stream_output(),
                  input_id: make_stream_input()}

    #6.sh LEASH OOM with -I flag 
    #TODO: tune it -m 
    

    if is_s3:
        flag_option_list = [OptionWithIO("-i", input_id),
                        OptionWithIO("-o", output_id),
                        Flag("-I"),
                        OptionWithIO("-m", ArgStringType(Arg.string_to_arg("1G"))),
                        OptionWithIO("-b", ArgStringType(Arg.string_to_arg(str(config.config['runtime']['dgsh_buffer_size']))))]
    else: #ec2 

        flag_option_list = [OptionWithIO("-i", input_id),
                        OptionWithIO("-o", output_id),
                        # Flag("-I"),
                        # OptionWithIO("-m", ArgStringType(Arg.string_to_arg("2G"))),
                        OptionWithIO("-b", ArgStringType(Arg.string_to_arg(str(config.config['runtime']['dgsh_buffer_size']))))]

        # Flag("-f"),
        # OptionWithIO("-T", ArgStringType(Arg.string_to_arg("/tmp"))),
        # OptionWithIO("-t", ArgStringType(Arg.string_to_arg(""))),

    # print(flag_option_list)

    cmd_inv_with_io_vars = CommandInvocationWithIOVars(
        cmd_name=dgsh_tee_bin,
        flag_option_list=flag_option_list,
        operand_list=[],
        implicit_use_of_streaming_input=None,
        implicit_use_of_streaming_output=None,
        access_map=access_map)
    return DGSHTee(cmd_inv_with_io_vars)
