from definitions.ir.dfg_node import *

class RemotePipe(DFGNode):
    def __init__(self, inputs, outputs, com_name, com_category,
                 com_options = [], com_redirs = [], com_assignments=[]):
        super().__init__(inputs, outputs, com_name, com_category,
                         com_options=com_options, 
                         com_redirs=com_redirs, 
                         com_assignments=com_assignments)

    def is_remote_read(self):
        com_name = self.com_name.opt_serialize()
        read_com = config.config['runtime']['remote_read_binary']
        return read_com in com_name

def make_remote_pipe(inputs, outputs, host_ip, port, is_remote_read, id):
    com_category = "pure"
    options = []
    opt_count = 0

    if is_remote_read:
        remote_pipe_bin = os.path.join(config.PASH_TOP, config.config['runtime']['remote_read_binary'])
    else:
        remote_pipe_bin = os.path.join(config.PASH_TOP, config.config['runtime']['remote_write_binary'])

    com_name = Arg(string_to_argument(remote_pipe_bin))

    options.append((opt_count, Arg(string_to_argument(f"--addr {host_ip}:{port}"))))
    options.append((opt_count + 1, Arg(string_to_argument(f"--id {id}"))))

    return RemotePipe(inputs,
               outputs,
               com_name, 
               com_category,
               com_options=options)
