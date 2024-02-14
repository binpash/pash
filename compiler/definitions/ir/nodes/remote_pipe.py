from definitions.ir.dfg_node import *
from pash_annotations.datatypes.CommandInvocationWithIOVars import CommandInvocationWithIOVars
from pash_annotations.datatypes.BasicDatatypes import ArgStringType, Flag
from pash_annotations.datatypes.BasicDatatypesWithIO import OptionWithIO
from pash_annotations.datatypes.AccessKind import make_stream_input, make_stream_output


class RemotePipe(DFGNode):
    def __init__(self, cmd_invocation_with_io_vars, addr=None, uuid=None):
        super().__init__(cmd_invocation_with_io_vars)
        self.addr = addr
        self.uuid = uuid
        

    def add_debug_flag(self):
        self.cmd_invocation_with_io_vars.flag_option_list.append(Flag("-d"))

    def is_remote_read(self):
        cmd_name = self.cmd_invocation_with_io_vars.cmd_name
        read_com = config.config['runtime']['remote_read_binary']
        return read_com in cmd_name
    
    def get_host(self):
        return self.addr[0]

    def get_local_host(self):
        for idx, option in enumerate(self.cmd_invocation_with_io_vars.flag_option_list):
            if "--localAddr" == option.option_name:
                return option.option_arg.name.opt_serialize().split(":")[0]
        return ""
    
    def get_uuid(self):
        return self.uuid

    def set_addr(self, host_ip, port):
        for idx, option in enumerate(self.cmd_invocation_with_io_vars.flag_option_list):
            if "--addr" == option.option_name:
                # Replace with new addr_option
                self.cmd_invocation_with_io_vars.flag_option_list[idx] = OptionWithIO("--addr", ArgStringType(Arg.string_to_arg(f"{host_ip}:{port}")))
                # Update metadata
                self.addr = (host_ip, port)

    def set_local_addr(self, host, port):
        for idx, option in enumerate(self.cmd_invocation_with_io_vars.flag_option_list):
            if "--localAddr" == option.option_name:
                # Replace with new addr_option
                self.cmd_invocation_with_io_vars.flag_option_list[idx] = OptionWithIO("--localAddr", ArgStringType(Arg.string_to_arg(f"{host}:{port}")))
            return
        self.cmd_invocation_with_io_vars.flag_option_list.append(OptionWithIO("--localAddr", ArgStringType(Arg.string_to_arg(f"{host}:{port}"))))

def make_remote_pipe(inputs, outputs, host_ip, port, is_remote_read, id):
    if is_remote_read:
        remote_pipe_bin = os.path.join(
            config.DISH_TOP, config.config['runtime']['remote_read_binary'])
    else:
        remote_pipe_bin = os.path.join(
            config.DISH_TOP, config.config['runtime']['remote_write_binary'])

    options = []
    if is_remote_read:
        options.append(OptionWithIO("--localAddr", ArgStringType(Arg.string_to_arg(f"{host_ip}:{port}"))))
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
        ),
        addr=(host_ip, port),
        uuid=id
    )
