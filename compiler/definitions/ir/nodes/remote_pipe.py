from definitions.ir.dfg_node import *
from pash_annotations.datatypes.CommandInvocationWithIOVars import CommandInvocationWithIOVars
from pash_annotations.datatypes.BasicDatatypes import ArgStringType, Flag
from pash_annotations.datatypes.BasicDatatypesWithIO import OptionWithIO
from pash_annotations.datatypes.AccessKind import make_stream_input, make_stream_output


class RemotePipe(DFGNode):
    def __init__(self, inputs, outputs, com_name, com_category,
                com_options=[], com_redirs=[], com_assignments=[], 
                addr=None, uuid=None):
        super().__init__(inputs, outputs, com_name, com_category,
                        com_options=com_options,
                        com_redirs=com_redirs,
                        com_assignments=com_assignments)
        self.addr = addr
        self.uuid = uuid
        

    def add_debug_flag(self):
        self.cmd_invocation_with_io_vars.flag_option_list.append(Flag("-d"))


    def is_remote_read(self):
        com_name = self.com_name.opt_serialize()
        read_com = config.config['runtime']['remote_read_binary']
        return read_com in com_name
    
    def get_host(self):
        return self.addr[0]
    
    def get_uuid(self):
        return self.uuid

    def set_addr(self, host_ip, port):
        for idx, option in enumerate(self.com_options):
            if "--addr" in option[1].opt_serialize():
                # Replace with new addr_option
                self.com_options[idx] = (idx, Arg(string_to_argument(f"--addr {host_ip}:{port}")))
                # Update metadata
                self.addr = (host_ip, port)

    def set_addr_conditional(self, uuid, host_ip, port):
        # when the remote_pipe has id specified,
        # update with new host_ip, port
        for idx, option in enumerate(self.com_options):
            if "--addr" in option[1].opt_serialize():
                addr = option[1].opt_serialize().split(' ')[1]
                if self.get_uuid() == uuid:
                    # Replace with new addr_option
                    self.com_options[idx] = (idx, Arg(string_to_argument(f"--addr {host_ip}:{port}")))
                    # Update metadata
                    self.addr = (host_ip, port)


def make_remote_pipe(inputs, outputs, host_ip, port, is_remote_read, id):
    if is_remote_read:
        remote_pipe_bin = os.path.join(
            config.DISH_TOP, config.config['runtime']['remote_read_binary'])
    else:
        remote_pipe_bin = os.path.join(
            config.DISH_TOP, config.config['runtime']['remote_write_binary'])

    com_name = Arg(string_to_argument(remote_pipe_bin))

    options.append(
        (opt_count, Arg(string_to_argument(f"--addr {host_ip}:{port}"))))
    options.append((opt_count + 1, Arg(string_to_argument(f"--id {id}"))))

    return RemotePipe(inputs,
                      outputs,
                      com_name,
                      com_category,
                      com_options=options,
                      addr=(host_ip, port),
                      uuid=id
                      )
