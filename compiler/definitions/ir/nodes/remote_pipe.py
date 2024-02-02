from definitions.ir.dfg_node import *
from pash_annotations.datatypes.CommandInvocationWithIOVars import CommandInvocationWithIOVars
from pash_annotations.datatypes.BasicDatatypes import ArgStringType, Flag
from pash_annotations.datatypes.BasicDatatypesWithIO import OptionWithIO
from pash_annotations.datatypes.AccessKind import make_stream_input, make_stream_output


class RemotePipe(DFGNode):
    def __init__(self, cmd_invocation_with_io_vars):
        super().__init__(cmd_invocation_with_io_vars)

    def add_debug_flag(self):
        self.cmd_invocation_with_io_vars.flag_option_list.append(Flag("-d"))

    def add_kill_flag(self, addr):
        self.cmd_invocation_with_io_vars.flag_option_list.append(
            OptionWithIO("--kill", ArgStringType(Arg.string_to_arg(addr))))

    def is_remote_read(self):
        cmd_name = self.cmd_invocation_with_io_vars.cmd_name
        read_com = config.config['runtime']['remote_read_binary']
        return read_com in cmd_name


def make_remote_pipe(inputs, outputs, host_ip, port, is_remote_read, id):
    if is_remote_read:
        remote_pipe_bin = os.path.join(
            config.DISH_TOP, config.config['runtime']['remote_read_binary'])
    else:
        remote_pipe_bin = os.path.join(
            config.DISH_TOP, config.config['runtime']['remote_write_binary'])

    options = []
    options.append(OptionWithIO("--addr", ArgStringType(Arg.string_to_arg(f"{host_ip}:{port}"))))
    options.append(OptionWithIO("--id", ArgStringType(Arg.string_to_arg(str(id)))))

    access_map = access_map = {
        **{input_id: make_stream_input() for input_id in inputs},
        **{output_id: make_stream_output() for output_id in outputs}
    }

    return RemotePipe(
        CommandInvocationWithIOVars(
            cmd_name=remote_pipe_bin,
            flag_option_list=options,
            operand_list=[],
            # This would only work if there is at most single input and output
            implicit_use_of_streaming_input = inputs[0] if inputs else None,
            implicit_use_of_streaming_output = outputs[0] if outputs else None,
            access_map=access_map
        )
    )
