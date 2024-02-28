from definitions.ir.dfg_node import *

class RemotePipe(DFGNode):
    def __init__(self, inputs, outputs, com_name, com_category,
                 com_options=[], com_redirs=[], com_assignments=[],
                 addr=None, uuid=None):
        super().__init__(inputs, outputs, com_name, com_category,
                         com_options=com_options,
                         com_redirs=com_redirs,
                         com_assignments=com_assignments
                        )
        self.addr = addr
        self.uuid = uuid


    def add_debug_flag(self):
        opt_count = len(self.com_options)
        self.com_options.append((opt_count, Arg(string_to_argument(f"-d"))))

    def is_remote_read(self):
        com_name = self.com_name.opt_serialize()
        read_com = config.config['runtime']['remote_read_binary']
        return read_com in com_name

    def get_host(self):
        return self.addr[0]

    def get_local_host(self):
        for idx, option in enumerate(self.com_options):
            if "--localAddr" in option[1].opt_serialize():
                return option[1].opt_serialize().split(' ')[1].split(':')[0]
        return ""

    def get_uuid(self):
        return self.uuid

    def set_addr(self, host, port):
        for idx, option in enumerate(self.com_options):
            if "--addr" in option[1].opt_serialize():
                # Replace with new addr_option
                self.com_options[idx] = (idx, Arg(string_to_argument(f"--addr {host}:{port}")))
                # Update metadata
                self.addr = (host, port)

    def set_local_addr(self, host, port):
        for idx, option in enumerate(self.com_options):
            if "--localAddr" in option[1].opt_serialize():
                # Replace with new addr_option
                self.com_options[idx] = (idx, Arg(string_to_argument(f"--localAddr {host}:{port}")))
            return
        self.com_options.append(Arg(string_to_argument(f"--localAddr {host}:{port}")))




def make_remote_pipe(inputs, outputs, host, port, is_remote_read, id):
    com_category = "pure"
    options = []

    if is_remote_read:
        remote_pipe_bin = os.path.join(
            config.DISH_TOP, config.config['runtime']['remote_read_binary'])
    else:
        remote_pipe_bin = os.path.join(
            config.DISH_TOP, config.config['runtime']['remote_write_binary'])

    com_name = Arg(string_to_argument(remote_pipe_bin))

    if is_remote_read:
        options.append((len(options), Arg(string_to_argument(f"--localAddr {host}:{port}"))))
    options.append(
        (len(options), Arg(string_to_argument(f"--addr {host}:{port}"))))
    options.append((len(options), Arg(string_to_argument(f"--id {id}"))))

    return RemotePipe(inputs,
                      outputs,
                      com_name,
                      com_category,
                      com_options=options,
                      addr=(host, port),
                      uuid=id)
