from definitions.ir.dfg_node import *


class RemotePipe(DFGNode):
    def __init__(self, inputs, outputs, com_name, com_category,
                 com_options=[], com_redirs=[], com_assignments=[]):
        super().__init__(inputs, outputs, com_name, com_category,
                         com_options=com_options,
                         com_redirs=com_redirs,
                         com_assignments=com_assignments)

    def add_debug_flag(self):
        opt_count = len(self.com_options)
        self.com_options.append((opt_count, Arg(string_to_argument(f"-d"))))

    def add_kill_flag(self, addr):
        opt_count = len(self.com_options)
        self.com_options.append((opt_count, Arg(string_to_argument(f"--kill {addr}"))))

    def add_manager_server_addr(self, host, port):
        opt_count = len(self.com_options)
        self.com_options.append((opt_count, Arg(string_to_argument(f"--managerAddr {host}:{port}"))))

    def is_remote_read(self):
        com_name = self.com_name.opt_serialize()
        read_com = config.config['runtime']['remote_read_binary']
        return read_com in com_name

    def get_manager_server_host(self):
        for idx, option in enumerate(self.com_options):
            if "--managerAddr" in option[1].opt_serialize():
                return option[1].opt_serialize().split(' ')[1].split(':')[0]
        return ""


def make_remote_pipe(inputs, outputs, host_ip, port, is_remote_read, id):
    com_category = "pure"
    options = []
    opt_count = 0

    if is_remote_read:
        remote_pipe_bin = os.path.join(
            config.DISH_TOP, config.config['runtime']['remote_read_binary'])
    else:
        remote_pipe_bin = os.path.join(
            config.DISH_TOP, config.config['runtime']['remote_write_binary'])

    com_name = Arg(string_to_argument(remote_pipe_bin))
    # serverAddr is not necessary - just passing host_ip would be enough
    options.append(
        (opt_count, Arg(string_to_argument(f"--addr {host_ip}:{port}"))))
    options.append((opt_count + 1, Arg(string_to_argument(f"--id {id}"))))

    return RemotePipe(inputs,
                      outputs,
                      com_name,
                      com_category,
                      com_options=options)
