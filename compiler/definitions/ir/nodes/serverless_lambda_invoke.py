from definitions.ir.dfg_node import *
from pash_annotations.datatypes.CommandInvocationWithIOVars import CommandInvocationWithIOVars
# from pash_annotations.datatypes.AccessKind import make_stream_input, make_stream_output
from pash_annotations.datatypes.BasicDatatypes import Operand

class ServerlessLambdaInvoke(DFGNode):
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
    
    def change_instance(self, instance):
        self.cmd_invocation_with_io_vars.operand_list[-1] = Operand(Arg.string_to_arg(str(instance)))

def make_serverless_lambda_invoke(key):
    lambda_invoke_bin = os.path.join('aws','invoke-lambda.py')
    cmd_inv_with_io_vars = CommandInvocationWithIOVars(
        cmd_name="python3 "+lambda_invoke_bin,
        flag_option_list=[],
        operand_list=[Operand(Arg.string_to_arg(str(key))), Operand(Arg.string_to_arg("$1")), Operand(Arg.string_to_arg("lambda"))],
        implicit_use_of_streaming_input=None,
        implicit_use_of_streaming_output=None,
        access_map={})
    return ServerlessLambdaInvoke(cmd_inv_with_io_vars)
